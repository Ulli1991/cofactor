"""Tests of lib/extract.py: every failure rule, the manifest, variants. Run: python3 -m unittest discover tests"""
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
import extract  # noqa: E402

SRC = """\
#include <math.h>
/* start of constants */
static const double A = 2.0;
static const double B = 3.0;
/* end of constants */

double f(double x)
{
    const char *s = "a } brace in a string";
    /* a } brace in a comment */
    if(x > 0) {return A * x;}
    return B * x;
}

static inline int g(int n) {
    return n + 1;
}

#ifdef FEATURE
double h(double x)
{
    return x;
}
#endif
"""


class Base(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        (self.root / "src").mkdir()
        (self.root / "src" / "code.c").write_text(SRC)
        self.X = extract.Extractor(self.root / "src", self.root / "gen", self.root / "man")

    def tearDown(self):
        self.tmp.cleanup()

    def out(self, name):
        return (self.root / "gen" / name).read_text()


class TestBlock(Base):
    def test_block_ok(self):
        self.X.block("c.h", "code.c", "/* start of constants */", "/* end of constants */")
        t = self.out("c.h")
        self.assertIn("A = 2.0", t)
        self.assertNotIn("end of constants", t)

    def test_include_stop(self):
        self.X.block("c.h", "code.c", "/* start of constants */", "/* end of constants */", include_stop=True)
        self.assertIn("end of constants", self.out("c.h"))

    def test_missing_start_marker(self):
        with self.assertRaises(extract.ExtractError):
            self.X.block("c.h", "code.c", "no such marker", "/* end")

    def test_duplicate_start_marker(self):
        with self.assertRaises(extract.ExtractError):
            self.X.block("c.h", "code.c", "static const double", "/* end")

    def test_missing_stop_marker(self):
        with self.assertRaises(extract.ExtractError):
            self.X.block("c.h", "code.c", "/* start of constants */", "never appears")

    def test_missing_file(self):
        with self.assertRaises(extract.ExtractError):
            self.X.block("c.h", "nope.c", "a", "b")

    def test_start_after(self):
        self.X.block("r.h", "code.c", "    return n", "}", start_after="static inline int g(")
        self.assertIn("n + 1", self.out("r.h"))

    def test_start_after_still_needs_a_unique_marker(self):
        # "    return" occurs in g and in h after the start_after marker: ambiguous, so it must refuse
        with self.assertRaises(extract.ExtractError):
            self.X.block("r.h", "code.c", "    return", "}", start_after="static inline int g(")

    def test_preprocessor_repair_declared(self):
        self.X.block("h.h", "code.c", "#ifdef FEATURE", "    return x;", expect_closed=1)
        self.assertTrue(self.out("h.h").rstrip().endswith("#endif"))

    def test_preprocessor_repair_undeclared(self):
        with self.assertRaises(extract.ExtractError):
            self.X.block("h.h", "code.c", "#ifdef FEATURE", "    return x;")

    def test_cut_last_and_wrap(self):
        self.X.block("c.h", "code.c", "/* start of constants */", "/* end of constants */", cut_last="static const double B",
                     prefix="// pre\n", suffix="// post\n")
        t = self.out("c.h")
        self.assertTrue(t.startswith("// pre\n") and t.endswith("// post\n"))
        self.assertNotIn("B = 3.0", t)


class TestFunctions(Base):
    def test_function_braces_in_strings_and_comments(self):
        self.X.functions("f.h", [("code.c", r"^double f\(", "double f(double x)")])
        t = self.out("f.h")
        self.assertIn("return B * x;", t)
        self.assertTrue(t.rstrip().endswith("}"))
        self.assertNotIn("static inline int g", t)

    def test_function_brace_on_signature_line(self):
        self.X.functions("g.h", [("code.c", r"^static inline int g\(", "static int g2(int {arg})")])
        t = self.out("g.h")
        self.assertIn("static int g2(int n) {", t)
        self.assertIn("return n + 1;", t)

    def test_function_not_found(self):
        with self.assertRaises(extract.ExtractError):
            self.X.functions("x.h", [("code.c", r"^double nothere\(", "double nothere(void)")])

    def test_literal_lines(self):
        self.X.functions("m.h", ["#ifdef X", ("code.c", r"^double f\(", "double f(double x)"), "#endif"])
        t = self.out("m.h")
        self.assertTrue(t.startswith("#ifdef X") and "#endif" in t)


class TestVariantAndManifest(Base):
    def test_variant_ok(self):
        self.X.functions("f.h", [("code.c", r"^double f\(", "double f(double x)")])
        self.X.variant("f.h", "f_fixed.h", [("return B * x;", "return 4.0 * x;")])
        self.assertIn("return 4.0 * x;", self.out("f_fixed.h"))
        self.assertIn("return B * x;", self.out("f.h"))

    def test_variant_old_text_must_occur_once(self):
        self.X.functions("f.h", [("code.c", r"^double f\(", "double f(double x)")])
        with self.assertRaises(extract.ExtractError):
            self.X.variant("f.h", "v.h", [("x", "y")])
        with self.assertRaises(extract.ExtractError):
            self.X.variant("f.h", "v.h", [("not in the block", "y")])

    def test_manifest_created_unchanged_changed_accept(self):
        def run(accept=False):
            X = extract.Extractor(self.root / "src", self.root / "gen", self.root / "man")
            X.block("c.h", "code.c", "/* start of constants */", "/* end of constants */")
            return X.finish(accept=accept)
        self.assertEqual(run(), 0)  # created
        self.assertEqual(run(), 0)  # unchanged
        p = self.root / "src" / "code.c"
        p.write_text(p.read_text().replace("A = 2.0", "A = 2.5"))
        self.assertEqual(run(), 3)  # target changed
        self.assertEqual(run(accept=True), 0)
        self.assertEqual(run(), 0)  # new manifest holds


if __name__ == "__main__":
    unittest.main()
