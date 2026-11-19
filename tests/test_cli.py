"""Tests of lib/cofactor_cli.py pieces that need no compiler: placeholders, value inheritance, project discovery,
the batch guard. Run: python3 -m unittest discover tests"""
import os
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
import cofactor_cli as cli  # noqa: E402


class TestSubst(unittest.TestCase):
    def test_known_unknown_and_shell(self):
        ctx = {"a": "x", "b": "{a}/y", "words": "P Q=(1/3)"}
        self.assertEqual(cli.subst("{b}", ctx), "x/y")
        self.assertEqual(cli.subst("${HOME} {a}", ctx), "${HOME} x")
        self.assertEqual(cli.subst("{words:q}", ctx), "P 'Q=(1/3)'")
        with self.assertRaises(cli.Fail):
            cli.subst("{nope}", ctx)


class TestConfigDefines(unittest.TestCase):
    def test_config_sh_lines(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "Config.sh"
            p.write_text("MAGNETIC\n# comment line\nEOS_GAMMA=(5.0/3.0)   % trailing comment\n\nKERNEL_FUNCTION=9 # x\n  HYDRO_MESHLESS_FINITE_MASS\n")
            self.assertEqual(cli.config_defines(p), ["MAGNETIC", "EOS_GAMMA=(5.0/3.0)", "KERNEL_FUNCTION=9", "HYDRO_MESHLESS_FINITE_MASS"])
            with self.assertRaises(cli.Fail):
                cli.config_defines(Path(d) / "missing.sh")


class Project(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        r = Path(self.tmp.name)
        self.proj = r / "proj"
        c = self.proj / ".cofactor" / "checks" / "chk"
        c.mkdir(parents=True)
        (self.proj / ".cofactor" / "project.ini").write_text("[project]\nname = proj\n[defaults]\nstd = c++17\nharness_args = {cfg}\n[vars]\nv = project\n")
        (c / "target.ini").write_text("[target]\ndescription = d\n[vars]\nw = check\n[config:a]\ndefines = X=1\n[config:b]\nstd = c++20\nexpect = FAIL\n")
        (self.proj / "sub").mkdir()
        self.site = r / "site.ini"
        self.site.write_text("[toolchain]\nenzyme_plugin = /none\n[vars]\nv = site\n[batch]\nrequire_batch = yes\njob_env = COFACTOR_TEST_JOB\n[paths]\nwork = %s\n" % (r / "work"))
        self.old = (cli.START, os.environ.get("COFACTOR_SITE"), os.environ.pop("COFACTOR_TEST_JOB", None))
        os.environ["COFACTOR_SITE"] = str(self.site)
        cli.START = self.proj / "sub"

    def tearDown(self):
        cli.START = self.old[0]
        if self.old[1] is None:
            os.environ.pop("COFACTOR_SITE", None)
        else:
            os.environ["COFACTOR_SITE"] = self.old[1]
        self.tmp.cleanup()

    def test_project_found_from_subdirectory(self):
        self.assertEqual(cli.find_project(cli.START), self.proj.resolve())

    def test_target_and_inheritance(self):
        T = cli.load_target("chk")
        self.assertEqual(T.name, "proj/chk")
        self.assertEqual(cli.configs(T), ["a", "b"])
        self.assertEqual(cli.tval(T, "a", "std"), "c++17")   # project default
        self.assertEqual(cli.tval(T, "b", "std"), "c++20")   # config wins
        self.assertEqual(cli.tval(T, "b", "expect", "PASS"), "FAIL")
        S = cli.load_site()
        ctx = cli.context(S, T, "a", self.proj)
        self.assertEqual(ctx["v"], "site")                    # site vars override project vars
        self.assertEqual(ctx["w"], "check")
        self.assertEqual(ctx["defines"], "X=1")               # config keys are placeholders
        self.assertEqual(Path(ctx["project"]), self.proj.resolve())

    def test_default_src_is_project_root(self):
        T = cli.load_target("chk")
        self.assertEqual(cli.resolve_src(T, "a", None, True), self.proj.resolve())

    def test_batch_guard(self):
        S = cli.load_site()
        with self.assertRaises(cli.Fail):
            cli.require_batch(S, "build")
        os.environ["COFACTOR_TEST_JOB"] = "1"
        try:
            cli.require_batch(S, "build")
        finally:
            os.environ.pop("COFACTOR_TEST_JOB")

    def test_workdir_tagged_by_source(self):
        S = cli.load_site()
        T = cli.load_target("chk")
        a = cli.workdir(S, T, "a", self.proj / "x")
        b = cli.workdir(S, T, "a", self.proj / "y")
        self.assertNotEqual(a, b)


if __name__ == "__main__":
    unittest.main()
