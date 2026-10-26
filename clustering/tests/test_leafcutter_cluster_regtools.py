"""Unit tests for clustering/leafcutter_cluster_regtools.py.

Run from the repository root with:

    python3 -B -m unittest discover -s clustering/tests

Inputs are small regtools .junc files written by each test. Expected outputs are worked
out by hand from the clustering rules in the script: an intron is read from a junction as
(chromStart + blockSizes[0], chromEnd - blockSizes[1] + 1), pooled across samples per
chromosome and strand, kept when it has at least 3 pooled reads, grouped by overlap and
then by shared splice sites, and filtered by -m and -p.
"""

import contextlib
import gzip
import importlib.util
import io
import os
import shutil
import subprocess
import sys
import tempfile
import types
import unittest
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.join(os.path.dirname(HERE), "leafcutter_cluster_regtools.py")


def load_script():
    spec = importlib.util.spec_from_file_location("leafcutter_cluster_regtools", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    # keep __pycache__ out of clustering/
    dont_write = sys.dont_write_bytecode
    sys.dont_write_bytecode = True
    try:
        spec.loader.exec_module(module)
    finally:
        sys.dont_write_bytecode = dont_write
    return module


lc = load_script()


def junc(chrom, start, end, count, strand="+"):
    """Return a 12-column regtools line that the script reads as the intron start:end.

    With 10 bp blocks on both sides the BED interval is [start - 10, end + 9).
    """
    bed_start, bed_end = start - 10, end + 9
    return "%s\t%d\t%d\tJUNC\t%d\t%s\t%d\t%d\t255,0,0\t2\t10,10\t0,%d\n" % (
        chrom, bed_start, bed_end, count, strand, bed_start, bed_end,
        bed_end - bed_start - 10)


MULTI_BLOCK = "chr1\t90\t420\tJUNC\t50\t+\t90\t420\t255,0,0\t3\t10,10,10\t0,150,320\n"
SHORT = "chr1\t100\t200\n"


def write_text(path, text):
    data = text.encode("utf-8")
    if path.endswith(".gz"):
        with gzip.open(path, "wb") as f:
            f.write(data)
    else:
        with open(path, "wb") as f:
            f.write(data)


def read_text(path):
    if path.endswith(".gz"):
        with gzip.open(path, "rb") as f:
            return f.read().decode("utf-8")
    with open(path, "rb") as f:
        return f.read().decode("utf-8")


# Three samples (plain, gzipped, plain). s1 also carries lines that must not change any
# count: a nonstandard chromosome, a "?" strand, an intron longer than the default -l, a
# three-block line and a short line.
MAIN_SAMPLES = [
    ("s1.junc", [
        junc("chr1", 100, 201, 12),
        junc("chr1", 100, 301, 3),
        junc("chr1", 250, 301, 9),
        junc("chrM", 100, 201, 50),
        junc("chrM", 100, 301, 50),
        junc("chr1", 100, 201, 50, "?"),
        junc("chr1", 100, 200201, 50),
        MULTI_BLOCK,
        SHORT,
        junc("chr1", 950, 1101, 10),
        junc("chr1", 1000, 1101, 4),
        junc("chr1", 2000, 2101, 1),
        junc("chr1", 2000, 2201, 5),
        junc("chr1", 500, 601, 6, "-"),
        junc("chr1", 500, 701, 6, "-"),
        junc("chr2", 300, 401, 2),
        junc("chr2", 300, 501, 1),
    ]),
    ("s2.junc.gz", [
        junc("chr1", 500, 601, 5, "-"),
        junc("chr1", 100, 201, 1),
        junc("chr1", 100, 301, 7),
        junc("chr1", 250, 301, 4),
        junc("chr1", 1000, 1101, 8),
        junc("chr1", 2000, 2101, 1),
        junc("chr1", 500, 701, 4, "-"),
        junc("chr2", 300, 401, 1),
        junc("chr2", 300, 501, 3),
    ]),
    ("s3.junc", [
        junc("chr1", 100, 201, 2),
        junc("chr1", 950, 1101, 1),
    ]),
]

# Pooled reads: 100:201 = 12+1+2, 100:301 = 3+7, 250:301 = 9+4, 950:1101 = 10+1,
# 1000:1101 = 4+8, 2000:2101 = 2 (< 3, dropped, which leaves 2000:2201 alone and unwritten),
# 500:601 = 6+5, 500:701 = 6+4, 300:401 = 3, 300:501 = 4. With -m 10 the chr2 cluster
# (7 reads) is dropped. Rows of the per-sample files are ordered by the string form of the
# coordinates within a cluster, so 1000:1101 comes before 950:1101.
MAIN_POOLED = (
    "chr1:+ 100:201:15 100:301:10 250:301:13 \n"
    "chr1:+ 950:1101:11 1000:1101:12 \n"
    "chr1:- 500:601:11 500:701:10 \n"
    "chr2:+ 300:401:3 300:501:4 \n"
)
MAIN_REFINED = (
    "chr1:+ 100:201:15 100:301:10 250:301:13 \n"
    "chr1:+ 950:1101:11 1000:1101:12 \n"
    "chr1:- 500:601:11 500:701:10 \n"
)
MAIN_SORTED = {
    "s1.junc": (
        "chrom s1\n"
        "chr1:100:201:clu_1_+ 12/24\n"
        "chr1:100:301:clu_1_+ 3/24\n"
        "chr1:250:301:clu_1_+ 9/24\n"
        "chr1:1000:1101:clu_2_+ 4/14\n"
        "chr1:950:1101:clu_2_+ 10/14\n"
        "chr1:500:601:clu_3_- 6/12\n"
        "chr1:500:701:clu_3_- 6/12\n"
    ),
    "s2.junc.gz": (
        "chrom s2\n"
        "chr1:100:201:clu_1_+ 1/12\n"
        "chr1:100:301:clu_1_+ 7/12\n"
        "chr1:250:301:clu_1_+ 4/12\n"
        "chr1:1000:1101:clu_2_+ 8/8\n"
        "chr1:950:1101:clu_2_+ 0/8\n"
        "chr1:500:601:clu_3_- 5/9\n"
        "chr1:500:701:clu_3_- 4/9\n"
    ),
    "s3.junc": (
        "chrom s3\n"
        "chr1:100:201:clu_1_+ 2/2\n"
        "chr1:100:301:clu_1_+ 0/2\n"
        "chr1:250:301:clu_1_+ 0/2\n"
        "chr1:1000:1101:clu_2_+ 0/1\n"
        "chr1:950:1101:clu_2_+ 1/1\n"
        "chr1:500:601:clu_3_- 0/0\n"
        "chr1:500:701:clu_3_- 0/0\n"
    ),
}
MAIN_COUNTS = (
    "chrom s1 s2 s3\n"
    "chr1:100:201:clu_1_+ 12/24 1/12 2/2\n"
    "chr1:100:301:clu_1_+ 3/24 7/12 0/2\n"
    "chr1:250:301:clu_1_+ 9/24 4/12 0/2\n"
    "chr1:1000:1101:clu_2_+ 4/14 8/8 0/1\n"
    "chr1:950:1101:clu_2_+ 10/14 0/8 1/1\n"
    "chr1:500:601:clu_3_- 6/12 5/9 0/0\n"
    "chr1:500:701:clu_3_- 6/12 4/9 0/0\n"
)
MAIN_NUMERS = (
    "s1 s2 s3\n"
    "chr1:100:201:clu_1_+ 12 1 2\n"
    "chr1:100:301:clu_1_+ 3 7 0\n"
    "chr1:250:301:clu_1_+ 9 4 0\n"
    "chr1:1000:1101:clu_2_+ 4 8 0\n"
    "chr1:950:1101:clu_2_+ 10 0 1\n"
    "chr1:500:601:clu_3_- 6 5 0\n"
    "chr1:500:701:clu_3_- 6 4 0\n"
)


def main_expected(rundir, prefix="main", with_pooling=True):
    """Return {file name in rundir: content} for a run on MAIN_SAMPLES with -m 10."""
    expected = {}
    if with_pooling:
        expected[prefix + "_pooled"] = MAIN_POOLED
        expected[prefix + "_refined"] = MAIN_REFINED
    sorted_names = ["%s.%s.sorted.gz" % (name, prefix) for name, _ in MAIN_SAMPLES]
    expected[prefix + "_sortedlibs"] = "".join(
        "%s/%s\n" % (rundir, name) for name in sorted_names)
    for (name, _), sorted_name in zip(MAIN_SAMPLES, sorted_names):
        expected[sorted_name] = MAIN_SORTED[name]
    expected[prefix + "_perind.counts.gz"] = MAIN_COUNTS
    expected[prefix + "_perind_numers.counts.gz"] = MAIN_NUMERS
    return expected


class TempDirTestCase(unittest.TestCase):

    def setUp(self):
        self.root = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.root)
        self.indir = os.path.join(self.root, "in")
        self.systmp = os.path.join(self.root, "systmp")
        os.mkdir(self.indir)
        os.mkdir(self.systmp)
        self.allowed = [self.indir]
        self.rundir = self.make_rundir("out")

    def make_rundir(self, name):
        path = os.path.join(self.root, name)
        os.mkdir(path)
        self.allowed.append(path)
        return path

    def write_juncs(self, samples):
        """Write [(file name, lines)] into the input directory and return the paths."""
        paths = []
        for name, lines in samples:
            path = os.path.join(self.indir, name)
            write_text(path, "".join(lines))
            paths.append(path)
        return paths

    def options(self, **kwargs):
        values = dict(juncfiles=None, outprefix="t", verbose=False, rundir=self.rundir,
                      maxintronlen=100000, minclureads=30, mincluratio=0.001,
                      cluster=None, nochromcheck=False, const=False)
        values.update(kwargs)
        return types.SimpleNamespace(**values)

    def call(self, func, *args, tmp_files_allowed=False):
        """Call func from self.root with the system temp dir moved under it.

        Checks that nothing lands outside the input and run directories and returns
        (stdout, stderr).
        """
        out, err = io.StringIO(), io.StringIO()
        cwd = os.getcwd()
        os.chdir(self.root)
        try:
            with mock.patch.object(tempfile, "tempdir", self.systmp), \
                    contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
                func(*args)
        finally:
            os.chdir(cwd)
        self.assert_no_stray_files(tmp_files_allowed)
        return out.getvalue(), err.getvalue()

    def output(self, name, rundir=None):
        return read_text(os.path.join(rundir or self.rundir, name))

    def assert_outputs(self, expected, rundir=None):
        rundir = rundir or self.rundir
        self.assertEqual(sorted(os.listdir(rundir)), sorted(expected))
        for name, content in expected.items():
            with self.subTest(file=name):
                self.assertEqual(self.output(name, rundir), content)

    def assert_no_stray_files(self, tmp_files_allowed=False):
        """Files sit in the input or a run directory; other directories are in the temp dir.

        self.root is also the working directory. The script leaves one empty directory in
        the temp dir per merge, and with more than 100 samples merged files in them too
        (tmp_files_allowed).
        """
        for dirpath, dirnames, filenames in os.walk(self.root):
            for name in dirnames:
                path = os.path.join(dirpath, name)
                if path != self.systmp and path not in self.allowed:
                    self.assertEqual(dirpath, self.systmp, "stray directory " + path)
            for name in filenames:
                path = os.path.join(dirpath, name)
                if not (tmp_files_allowed and path.startswith(self.systmp + os.sep)):
                    self.assertIn(dirpath, self.allowed, "stray file " + path)


class TestOverlaps(unittest.TestCase):

    def test_disjoint(self):
        self.assertFalse(lc.overlaps((1, 5), (6, 10)))
        self.assertFalse(lc.overlaps((6, 10), (1, 5)))

    def test_shared_endpoint_counts_as_overlap(self):
        self.assertTrue(lc.overlaps((1, 5), (5, 10)))
        self.assertTrue(lc.overlaps((5, 10), (1, 5)))

    def test_partial_and_contained(self):
        self.assertTrue(lc.overlaps((1, 6), (4, 10)))
        self.assertTrue(lc.overlaps((1, 10), (3, 4)))
        self.assertTrue(lc.overlaps((3, 4), (1, 10)))


class TestClusterIntervals(unittest.TestCase):

    def test_single_interval(self):
        clusters, ordered = lc.cluster_intervals([(10, 20)])
        self.assertEqual(clusters, [[(10, 20)]])
        self.assertEqual(ordered, [(10, 20)])

    def test_sorts_then_groups(self):
        clusters, ordered = lc.cluster_intervals([(30, 40), (10, 20), (15, 25)])
        self.assertEqual(clusters, [[(10, 20), (15, 25)], [(30, 40)]])
        self.assertEqual(ordered, [(10, 20), (15, 25), (30, 40)])

    def test_touching_intervals_join_and_gap_splits(self):
        clusters, _ = lc.cluster_intervals([(1, 10), (10, 20), (21, 30)])
        self.assertEqual(clusters, [[(1, 10), (10, 20)], [(21, 30)]])

    def test_extent_uses_largest_end_seen(self):
        """(50, 60) misses (10, 20) but lies inside (1, 100), so all three cluster."""
        clusters, _ = lc.cluster_intervals([(1, 100), (10, 20), (50, 60)])
        self.assertEqual(clusters, [[(1, 100), (10, 20), (50, 60)]])

    def test_chain_of_overlaps(self):
        clusters, _ = lc.cluster_intervals([(18, 30), (1, 10), (5, 20), (40, 50), (45, 60)])
        self.assertEqual(clusters, [[(1, 10), (5, 20), (18, 30)], [(40, 50), (45, 60)]])


class TestRefineLinked(unittest.TestCase):

    def test_shared_start_then_shared_end(self):
        introns = [((100, 200), 5), ((100, 300), 3), ((250, 300), 2)]
        self.assertEqual(lc.refine_linked(introns), [introns])

    def test_end_of_one_is_start_of_other(self):
        introns = [((100, 200), 1), ((200, 300), 1)]
        self.assertEqual(lc.refine_linked(introns), [introns])

    def test_order_follows_linking_passes(self):
        """(150, 300) only links through (200, 300), which joins first."""
        a, b, c = ((100, 200), 1), ((150, 300), 2), ((200, 300), 3)
        self.assertEqual(lc.refine_linked([a, b, c]), [[a, c, b]])

    def test_two_groups(self):
        a, b, c, d = ((100, 200), 1), ((100, 300), 2), ((400, 500), 3), ((400, 600), 4)
        self.assertEqual(lc.refine_linked([a, b, c, d]), [[a, b], [c, d]])

    def test_leading_unlinked_intron_is_own_group(self):
        s, a, b = ((50, 60), 1), ((100, 200), 2), ((100, 300), 3)
        self.assertEqual(lc.refine_linked([s, a, b]), [[s], [a, b]])


class TestRefineCluster(unittest.TestCase):

    def test_all_introns_pass(self):
        clu = [((100, 200), 10), ((100, 300), 10)]
        self.assertEqual(lc.refine_cluster(clu, 0.1, 1), [clu])

    def test_low_ratio_intron_removed(self):
        """5 of 100 reads is below 0.1; the other two still share position 100."""
        clu = [((100, 200), 50), ((100, 300), 45), ((150, 300), 5)]
        self.assertEqual(lc.refine_cluster(clu, 0.1, 1),
                         [[((100, 200), 50), ((100, 300), 45)]])

    def test_ratio_equal_to_cutoff_kept(self):
        clu = [((100, 200), 75), ((100, 300), 25)]
        self.assertEqual(lc.refine_cluster(clu, 0.25, 1), [clu])

    def test_read_cutoff_per_intron(self):
        clu = [((100, 200), 10), ((100, 300), 10), ((100, 400), 4)]
        self.assertEqual(lc.refine_cluster(clu, 0.0, 5),
                         [[((100, 200), 10), ((100, 300), 10)]])
        self.assertEqual(lc.refine_cluster(clu, 0.0, 4), [clu])

    def test_removing_bridge_splits_cluster(self):
        """Without 200:500 (4 of 124 reads, below 0.05) the rest form two separate pairs."""
        clu = [((100, 200), 30), ((100, 250), 30), ((200, 500), 4),
               ((400, 500), 30), ((400, 550), 30)]
        self.assertEqual(lc.refine_cluster(clu, 0.05, 1),
                         [[((100, 200), 30), ((100, 250), 30)],
                          [((400, 500), 30), ((400, 550), 30)]])
        self.assertEqual(lc.refine_cluster(clu, 0.03, 1), [clu])

    def test_overlapping_but_unlinked_groups_split(self):
        clu = [((100, 200), 10), ((100, 300), 10), ((150, 250), 10), ((150, 260), 10)]
        self.assertEqual(lc.refine_cluster(clu, 0.001, 1),
                         [[((100, 200), 10), ((100, 300), 10)],
                          [((150, 250), 10), ((150, 260), 10)]])

    def test_overlapping_unlinked_pair_gives_nothing(self):
        clu = [((100, 200), 10), ((150, 300), 10)]
        self.assertEqual(lc.refine_cluster(clu, 0.0, 1), [])

    def test_one_intron_left_gives_nothing(self):
        clu = [((100, 200), 1000), ((100, 300), 1)]
        self.assertEqual(lc.refine_cluster(clu, 0.01, 1), [])

    def test_no_intron_left_gives_nothing(self):
        clu = [((100, 200), 1), ((100, 300), 1)]
        self.assertEqual(lc.refine_cluster(clu, 0.0, 5), [])

    def test_output_order(self):
        """Input order does not matter; output follows linking order of the sorted introns."""
        a, b, c = ((100, 200), 10), ((150, 300), 10), ((200, 300), 10)
        self.assertEqual(lc.refine_cluster([c, b, a], 0.0, 1), [[a, c, b]])
        self.assertEqual(lc.refine_cluster([((100, 300), 5), ((100, 200), 5)], 0.0, 1),
                         [[((100, 200), 5), ((100, 300), 5)]])


class TestPoolJuncReads(TempDirTestCase):

    def pool(self, samples, **kwargs):
        paths = self.write_juncs(samples)
        out, err = self.call(lc.pool_junc_reads, paths, self.options(**kwargs))
        return self.output("t_pooled"), out, err, paths

    def test_intron_coordinates_from_block_sizes(self):
        """1000 + 7 = 1007; 1300 - 13 + 1 = 1288; 1400 - 20 + 1 = 1381."""
        lines = [
            "chr1\t1000\t1300\tJ1\t4\t+\t1000\t1300\t255,0,0\t2\t7,13\t0,287\n",
            "chr1\t1000\t1400\tJ2\t5\t+\t1000\t1400\t255,0,0\t2\t7,20\t0,380\n",
        ]
        pooled, _, _, _ = self.pool([("s1.junc", lines)])
        self.assertEqual(pooled, "chr1:+ 1007:1288:4 1007:1381:5 \n")

    def test_counts_summed_over_lines_and_samples(self):
        """100:201 has 1 + 1 + 1 = 3 reads (kept); 100:401 has 1 + 1 = 2 (dropped)."""
        pooled, _, _, _ = self.pool([
            ("s1.junc", [junc("chr1", 100, 201, 1), junc("chr1", 100, 301, 2),
                         junc("chr1", 100, 201, 1), junc("chr1", 100, 401, 1)]),
            ("s2.junc", [junc("chr1", 100, 201, 1), junc("chr1", 100, 301, 2),
                         junc("chr1", 100, 401, 1)]),
        ])
        self.assertEqual(pooled, "chr1:+ 100:201:3 100:301:4 \n")

    def test_row_order(self):
        """Chromosome/strand rows follow first appearance; introns are sorted numerically."""
        pooled, _, _, _ = self.pool([
            ("s1.junc", [junc("chr10", 5000, 5201, 4, "-"), junc("chr10", 5000, 5101, 3, "-"),
                         junc("chr2", 1000, 1101, 3), junc("chr2", 950, 1101, 3)]),
            ("s2.junc", [junc("chr1", 100, 301, 3), junc("chr1", 100, 201, 3),
                         junc("chr2", 950, 1101, 1), junc("chr10", 5000, 5101, 3),
                         junc("chr10", 5000, 5201, 3)]),
        ])
        self.assertEqual(pooled,
                         "chr10:- 5000:5101:3 5000:5201:4 \n"
                         "chr2:+ 950:1101:4 1000:1101:3 \n"
                         "chr1:+ 100:201:3 100:301:3 \n"
                         "chr10:+ 5000:5101:3 5000:5201:3 \n")

    def test_ambiguous_strand_skipped(self):
        pooled, _, _, _ = self.pool([("s1.junc", [
            junc("chr1", 100, 201, 5), junc("chr1", 100, 301, 5),
            junc("chr1", 100, 201, 50, "?"), junc("chr1", 100, 401, 50, "?"),
            junc("chr1", 100, 501, 50, "?")])])
        self.assertEqual(pooled, "chr1:+ 100:201:5 100:301:5 \n")

    def test_multi_block_lines_skipped_and_echoed(self):
        pooled, out, _, _ = self.pool([("s1.junc", [
            junc("chr1", 100, 201, 5), MULTI_BLOCK, junc("chr1", 100, 301, 5)])])
        self.assertEqual(pooled, "chr1:+ 100:201:5 100:301:5 \n")
        self.assertIn(MULTI_BLOCK.rstrip("\n"), out)

    def test_short_lines_skipped_with_message(self):
        pooled, _, err, paths = self.pool([("s1.junc", [
            junc("chr1", 100, 201, 5), SHORT, "\n", junc("chr1", 100, 301, 5)])])
        self.assertEqual(pooled, "chr1:+ 100:201:5 100:301:5 \n")
        self.assertEqual(err.count("Error in %s\n" % paths[0]), 2)

    CHROMS = ["chr1", "chrM", "chr22", "MT", "chrX", "chr23", "chrY", "23", "1",
              "chr1_KI270706v1_random", "22", "GL000220.1", "X", "chrx", "Y"]
    STANDARD = ["chr1", "chr22", "chrX", "chrY", "1", "22", "X", "Y"]

    def chrom_samples(self):
        lines = []
        for chrom in self.CHROMS:
            lines += [junc(chrom, 100, 201, 3), junc(chrom, 100, 301, 4)]
        return [("s1.junc", lines)]

    def test_chromosome_names_checked(self):
        pooled, _, _, _ = self.pool(self.chrom_samples())
        self.assertEqual(pooled, "".join("%s:+ 100:201:3 100:301:4 \n" % c
                                         for c in self.STANDARD))

    def test_nochromcheck_keeps_all_chromosomes(self):
        pooled, _, _, _ = self.pool(self.chrom_samples(), nochromcheck=True)
        self.assertEqual(pooled, "".join("%s:+ 100:201:3 100:301:4 \n" % c
                                         for c in self.CHROMS))

    def test_max_intron_length(self):
        """end - start of 100 is kept at -l 100 and 101 is not."""
        samples = [("s1.junc", [junc("chr1", 1000, 1050, 3), junc("chr1", 1000, 1100, 3),
                                junc("chr1", 1000, 1101, 3)])]
        pooled, _, _, _ = self.pool(samples, maxintronlen="100")
        self.assertEqual(pooled, "chr1:+ 1000:1050:3 1000:1100:3 \n")
        pooled, _, _, _ = self.pool(samples, maxintronlen="101")
        self.assertEqual(pooled, "chr1:+ 1000:1050:3 1000:1100:3 1000:1101:3 \n")

    def test_single_intron_clusters_not_written(self):
        """201 and 202 do not touch; chr2 has no intron with 3 reads."""
        pooled, _, _, _ = self.pool([("s1.junc", [
            junc("chr1", 100, 201, 5), junc("chr1", 202, 301, 5),
            junc("chr2", 100, 201, 2), junc("chr2", 100, 301, 2)])])
        self.assertEqual(pooled, "")

    def test_touching_introns_cluster(self):
        pooled, _, _, _ = self.pool([("s1.junc", [
            junc("chr1", 100, 201, 5), junc("chr1", 201, 301, 5)])])
        self.assertEqual(pooled, "chr1:+ 100:201:5 201:301:5 \n")

    def test_gzipped_input_matches_plain(self):
        lines = [junc("chr1", 100, 201, 4), junc("chr1", 100, 301, 5),
                 junc("chr1", 100, 401, 6, "-"), junc("chr1", 100, 501, 7, "-")]
        expected = "chr1:+ 100:201:4 100:301:5 \nchr1:- 100:401:6 100:501:7 \n"
        pooled, _, _, _ = self.pool([("p.junc", lines)])
        self.assertEqual(pooled, expected)
        pooled, _, _, _ = self.pool([("g.junc.gz", lines)])
        self.assertEqual(pooled, expected)


class TestRefineClusters(TempDirTestCase):

    def refine(self, pooled, **kwargs):
        write_text(os.path.join(self.rundir, "t_pooled"), pooled)
        self.call(lc.refine_clusters, self.options(**kwargs))
        return self.output("t_refined")

    def test_cluster_below_min_reads_dropped(self):
        refined = self.refine("chr1:+ 100:201:10 100:301:10 \n"
                              "chr1:+ 500:601:4 500:701:5 \n", minclureads="10")
        self.assertEqual(refined, "chr1:+ 100:201:10 100:301:10 \n")

    def test_low_ratio_intron_removed_and_cluster_split(self):
        pooled = "chr1:+ 100:200:30 100:250:30 200:500:4 400:500:30 400:550:30 \n"
        refined = self.refine(pooled, minclureads="3", mincluratio="0.05")
        self.assertEqual(refined, "chr1:+ 100:200:30 100:250:30 \n"
                                  "chr1:+ 400:500:30 400:550:30 \n")
        refined = self.refine(pooled, minclureads="3", mincluratio="0.03")
        self.assertEqual(refined, pooled)

    def test_unlinked_introns_form_separate_clusters(self):
        """120:280 shares no splice site with the others and is left out."""
        refined = self.refine(
            "chr1:+ 100:200:10 100:300:10 120:280:10 150:250:10 150:260:10 \n",
            minclureads="3")
        self.assertEqual(refined, "chr1:+ 100:200:10 100:300:10 \n"
                                  "chr1:+ 150:250:10 150:260:10 \n")

    def test_introns_written_in_linking_order(self):
        refined = self.refine("chr1:+ 100:200:10 150:300:10 200:300:10 \n",
                              minclureads="3")
        self.assertEqual(refined, "chr1:+ 100:200:10 200:300:10 150:300:10 \n")

    def test_cluster_with_one_surviving_intron_dropped(self):
        refined = self.refine("chr1:+ 100:201:1000 100:301:1 \n",
                              minclureads="1", mincluratio="0.01")
        self.assertEqual(refined, "")

    def test_const_keeps_single_intron_clusters(self):
        """With const, a one-intron line is kept when it has at least -m reads."""
        pooled = ("chr1:+ 100:201:40 \n"
                  "chr1:+ 300:401:10 \n"
                  "chr1:+ 500:601:10 500:701:10 \n"
                  "chr1:+ 900:1001:9 \n")
        refined = self.refine(pooled, minclureads="10", const=True)
        self.assertEqual(refined, "chr1:+ 100:201:40 \nchr1:+ 300:401:10 \n"
                                  "chr1:+ 500:601:10 500:701:10 \n")
        refined = self.refine(pooled, minclureads="10", const=False)
        self.assertEqual(refined, "chr1:+ 500:601:10 500:701:10 \n")


class TestSortJunctions(TempDirTestCase):

    def sort(self, samples, clusters, extra_paths=(), **kwargs):
        if kwargs.get("cluster") is None:
            write_text(os.path.join(self.rundir, "t_refined"), clusters)
        paths = self.write_juncs(samples) + list(extra_paths)
        return self.call(lc.sort_junctions, paths, self.options(**kwargs))

    def test_counts_and_totals_per_cluster(self):
        """Strands are kept apart, repeated lines add up, other introns are ignored."""
        self.sort([("s1.junc", [junc("chr1", 100, 201, 2), junc("chr1", 100, 301, 1),
                                junc("chr1", 100, 201, 7, "-"), junc("chr1", 100, 201, 3),
                                junc("chr1", 100, 401, 9), junc("chr1", 100, 301, 4, "?")])],
                  "chr1:+ 100:201:0 100:301:0 \nchr1:- 100:201:0 100:301:0 \n")
        self.assertEqual(self.output("s1.junc.t.sorted.gz"),
                         "chrom s1\n"
                         "chr1:100:201:clu_1_+ 5/6\n"
                         "chr1:100:301:clu_1_+ 1/6\n"
                         "chr1:100:201:clu_2_- 7/7\n"
                         "chr1:100:301:clu_2_- 0/7\n")
        self.assertEqual(self.output("t_sortedlibs"),
                         "%s/s1.junc.t.sorted.gz\n" % self.rundir)

    def test_cluster_absent_from_sample(self):
        self.sort([("s1.junc", [junc("chr1", 100, 201, 4)])],
                  "chr3:+ 100:201:0 100:301:0 \nchr1:+ 100:201:0 100:301:0 \n")
        self.assertEqual(self.output("s1.junc.t.sorted.gz"),
                         "chrom s1\n"
                         "chr3:100:201:clu_1_+ 0/0\n"
                         "chr3:100:301:clu_1_+ 0/0\n"
                         "chr1:100:201:clu_2_+ 4/4\n"
                         "chr1:100:301:clu_2_+ 0/4\n")

    def test_rows_sorted_as_strings_within_cluster(self):
        """Rows follow coordinate strings: "1000" < "950" < "99" and "1000" < "999"."""
        self.sort([("s1.junc", [junc("chr1", 950, 1101, 1), junc("chr1", 1000, 1101, 2),
                                junc("chr1", 99, 1101, 3), junc("chr2", 100, 999, 4),
                                junc("chr2", 100, 1000, 5)])],
                  "chr1:+ 950:1101:0 1000:1101:0 99:1101:0 \n"
                  "chr2:+ 100:999:0 100:1000:0 \n")
        self.assertEqual(self.output("s1.junc.t.sorted.gz"),
                         "chrom s1\n"
                         "chr1:1000:1101:clu_1_+ 2/6\n"
                         "chr1:950:1101:clu_1_+ 1/6\n"
                         "chr1:99:1101:clu_1_+ 3/6\n"
                         "chr2:100:1000:clu_2_+ 5/9\n"
                         "chr2:100:999:clu_2_+ 4/9\n")

    def test_cluster_file_option(self):
        """The cluster file is used instead of the run directory's t_refined."""
        clusters = os.path.join(self.indir, "clusters.txt")
        write_text(clusters, "chr1:+ 100:201:0 100:301:0 \n")
        write_text(os.path.join(self.rundir, "t_refined"), "chr1:+ 100:201:0 100:401:0 \n")
        self.sort([("s1.junc", [junc("chr1", 100, 201, 3), junc("chr1", 100, 301, 1),
                                junc("chr1", 100, 401, 7)])],
                  None, cluster=clusters)
        self.assertEqual(self.output("s1.junc.t.sorted.gz"),
                         "chrom s1\nchr1:100:201:clu_1_+ 3/4\nchr1:100:301:clu_1_+ 1/4\n")

    def test_missing_library_skipped(self):
        missing = os.path.join(self.indir, "missing.junc")
        self.sort([("s1.junc", [junc("chr1", 100, 201, 3)])],
                  "chr1:+ 100:201:0 100:301:0 \n", extra_paths=[missing])
        self.assertEqual(self.output("t_sortedlibs"),
                         "%s/s1.junc.t.sorted.gz\n" % self.rundir)
        self.assertEqual(sorted(os.listdir(self.rundir)),
                         ["s1.junc.t.sorted.gz", "t_refined", "t_sortedlibs"])

    def test_header_and_file_names(self):
        """The header keeps the file name up to the first ".junc"."""
        line = [junc("chr1", 100, 201, 1)]
        self.sort([("a.junc", line), ("b.junc.gz", line), ("c.txt", line)],
                  "chr1:+ 100:201:0 100:301:0 \n")
        names = ["a.junc.t.sorted.gz", "b.junc.gz.t.sorted.gz", "c.txt.t.sorted.gz"]
        self.assertEqual(self.output("t_sortedlibs"),
                         "".join("%s/%s\n" % (self.rundir, n) for n in names))
        for name, header in zip(names, ["a", "b", "c.txt"]):
            self.assertEqual(self.output(name),
                             "chrom %s\nchr1:100:201:clu_1_+ 1/1\n"
                             "chr1:100:301:clu_1_+ 0/1\n" % header)

    def test_chromosome_check(self):
        clusters = os.path.join(self.indir, "clusters.txt")
        write_text(clusters, "chrM:+ 100:201:0 100:301:0 \n")
        samples = [("s1.junc", [junc("chrM", 100, 201, 2), junc("chrM", 100, 301, 3)])]
        self.sort(samples, None, cluster=clusters)
        self.assertEqual(self.output("s1.junc.t.sorted.gz"),
                         "chrom s1\nchrM:100:201:clu_1_+ 0/0\nchrM:100:301:clu_1_+ 0/0\n")
        self.sort(samples, None, cluster=clusters, nochromcheck=True)
        self.assertEqual(self.output("s1.junc.t.sorted.gz"),
                         "chrom s1\nchrM:100:201:clu_1_+ 2/5\nchrM:100:301:clu_1_+ 3/5\n")

    def test_cluster_strand_matched_as_written(self):
        """A cluster file naming strand "?" collects the sample's "?" junctions.

        Pins the current behaviour: pooling never makes a "?" cluster.
        """
        clusters = os.path.join(self.indir, "clusters.txt")
        write_text(clusters, "chr1:? 100:201:0 100:301:0 \n")
        self.sort([("s1.junc", [junc("chr1", 100, 201, 2, "?"), junc("chr1", 100, 301, 3, "?"),
                                junc("chr1", 100, 201, 40, "+")])],
                  None, cluster=clusters)
        self.assertEqual(self.output("s1.junc.t.sorted.gz"),
                         "chrom s1\nchr1:100:201:clu_1_? 2/5\nchr1:100:301:clu_1_? 3/5\n")

    def test_short_and_multi_block_lines_ignored(self):
        out, err = self.sort([("s1.junc", [junc("chr1", 100, 201, 2), SHORT, MULTI_BLOCK,
                                           junc("chr1", 100, 301, 3)])],
                             "chr1:+ 100:201:0 100:301:0 \n")
        self.assertEqual(self.output("s1.junc.t.sorted.gz"),
                         "chrom s1\nchr1:100:201:clu_1_+ 2/5\nchr1:100:301:clu_1_+ 3/5\n")
        self.assertIn("Error in %s" % os.path.join(self.indir, "s1.junc"), err)
        self.assertIn(MULTI_BLOCK.rstrip("\n"), out)

    def test_repeated_path_read_once_per_listing(self):
        """A path listed twice gives one output file holding its counts twice.

        Pins the current behaviour.
        """
        path = self.write_juncs([("s1.junc", [junc("chr1", 100, 201, 2),
                                              junc("chr1", 100, 301, 3)])])[0]
        write_text(os.path.join(self.rundir, "t_refined"), "chr1:+ 100:201:0 100:301:0 \n")
        self.call(lc.sort_junctions, [path, path], self.options(verbose=True))
        self.assertEqual(self.output("t_sortedlibs"),
                         "%s/s1.junc.t.sorted.gz\n" % self.rundir)
        self.assertEqual(self.output("s1.junc.t.sorted.gz"),
                         "chrom s1\nchr1:100:201:clu_1_+ 4/10\nchr1:100:301:clu_1_+ 6/10\n")


class TestMergeFiles(TempDirTestCase):

    def test_columns_joined_line_by_line(self):
        paths = [os.path.join(self.indir, n) for n in ("a.txt", "b.gz", "c.txt")]
        write_text(paths[0], "chrom s1\nr1 1/2\nr2 1/2\n")
        write_text(paths[1], "chrom s2\nr1 0/3\nr2 3/3\n")
        write_text(paths[2], "chrom s3\nr1 4/4\nr2 0/4\n")
        fout = io.BytesIO()
        self.call(lc.merge_files, paths, fout, self.options(verbose=True))
        self.assertEqual(fout.getvalue().decode("utf-8"),
                         "chrom s1 s2 s3\nr1 1/2 0/3 4/4\nr2 1/2 3/3 0/4\n")


class TestMergeJunctions(TempDirTestCase):

    def merge(self, contents, tmp_files_allowed=False, **kwargs):
        names = []
        for i, text in enumerate(contents):
            name = os.path.join(self.rundir, "f%03d.sorted.gz" % i)
            write_text(name, text)
            names.append(name)
        write_text(os.path.join(self.rundir, "t_sortedlibs"), "".join(n + "\n" for n in names))
        self.call(lc.merge_junctions, self.options(**kwargs),
                  tmp_files_allowed=tmp_files_allowed)

    def test_two_files(self):
        self.merge(["chrom s1\nr1 1/2\nr2 1/2\n", "chrom s2\nr1 0/3\nr2 3/3\n"])
        self.assertEqual(self.output("t_perind.counts.gz"),
                         "chrom s1 s2\nr1 1/2 0/3\nr2 1/2 3/3\n")

    def test_single_file(self):
        self.merge(["chrom s1\nr1 1/2\nr2 1/2\n"])
        self.assertEqual(self.output("t_perind.counts.gz"), "chrom s1\nr1 1/2\nr2 1/2\n")

    def test_const_output_name(self):
        self.merge(["chrom s1\nr1 1/2\n", "chrom s2\nr1 2/2\n"], const=True)
        self.assertEqual(self.output("t_perind.constcounts.gz"), "chrom s1 s2\nr1 1/2 2/2\n")
        self.assertFalse(os.path.exists(os.path.join(self.rundir, "t_perind.counts.gz")))

    def test_more_than_one_hundred_files(self):
        """101 files are merged in two rounds; columns stay in list order."""
        n = 101
        self.merge(["chrom s%03d\nrow1 %d/1000\nrow2 1/%d\n" % (i, i, i + 1)
                    for i in range(n)], tmp_files_allowed=True)
        self.assertEqual(self.output("t_perind.counts.gz"),
                         "chrom " + " ".join("s%03d" % i for i in range(n)) + "\n"
                         "row1 " + " ".join("%d/1000" % i for i in range(n)) + "\n"
                         "row2 " + " ".join("1/%d" % (i + 1) for i in range(n)) + "\n")


class TestGetNumers(TempDirTestCase):

    def test_numerators(self):
        write_text(os.path.join(self.rundir, "t_perind.counts.gz"),
                   "chrom s1 s2\n"
                   "chr1:100:201:clu_1_+ 12/24 1/12\n"
                   "chr1:100:301:clu_1_+ 3/24 0/12\n")
        self.call(lc.get_numers, self.options())
        self.assertEqual(self.output("t_perind_numers.counts.gz"),
                         "s1 s2\n"
                         "chr1:100:201:clu_1_+ 12 1\n"
                         "chr1:100:301:clu_1_+ 3 0\n")

    def test_header_only(self):
        write_text(os.path.join(self.rundir, "t_perind.counts.gz"), "chrom s1 s2\n")
        self.call(lc.get_numers, self.options())
        self.assertEqual(self.output("t_perind_numers.counts.gz"), "s1 s2\n")

    def test_const_file_names(self):
        write_text(os.path.join(self.rundir, "t_perind.constcounts.gz"),
                   "chrom s1\nchr1:100:201:clu_1_- 7/9\n")
        self.call(lc.get_numers, self.options(const=True))
        self.assertEqual(self.output("t_perind_numers.constcounts.gz"),
                         "s1\nchr1:100:201:clu_1_- 7\n")


class TestMain(TempDirTestCase):

    def test_full_run(self):
        paths = self.write_juncs(MAIN_SAMPLES)
        self.call(lc.main, self.options(outprefix="main", minclureads=10, verbose=True),
                  paths)
        self.assert_outputs(main_expected(self.rundir))

    def test_cluster_file_skips_pooling(self):
        paths = self.write_juncs(MAIN_SAMPLES)
        clusters = os.path.join(self.indir, "clusters.txt")
        write_text(clusters, MAIN_REFINED)
        self.call(lc.main, self.options(outprefix="main", cluster=clusters), paths)
        self.assert_outputs(main_expected(self.rundir, with_pooling=False))


class TestCommandLine(TempDirTestCase):

    def run_script(self, *args, cwd=None):
        env = dict(os.environ)
        env["TMPDIR"] = self.systmp
        proc = subprocess.run([sys.executable, SCRIPT] + list(args), cwd=cwd or self.root,
                              env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                              universal_newlines=True)
        self.assert_no_stray_files()
        return proc

    def run_ok(self, *args, cwd=None):
        proc = self.run_script(*args, cwd=cwd)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        return proc

    def juncfile_list(self, samples):
        paths = self.write_juncs(samples)
        listing = os.path.join(self.indir, "juncfiles.txt")
        write_text(listing, "".join(p + "\n" for p in paths))
        return listing

    def test_main_scenario(self):
        listing = self.juncfile_list(MAIN_SAMPLES)
        s1 = os.path.join(self.indir, "s1.junc")
        styles = {
            "short": lambda r: ["-j", listing, "-r", r, "-o", "main", "-m", "10"],
            "long": lambda r: ["--juncfiles=" + listing, "--rundir=" + r,
                               "--outprefix=main", "--minclureads=10"],
        }
        for style, make_args in styles.items():
            with self.subTest(style=style):
                rundir = self.make_rundir("out-" + style)
                proc = self.run_ok(*make_args(rundir))
                self.assert_outputs(main_expected(rundir), rundir)
                self.assertIn(MULTI_BLOCK.rstrip("\n"), proc.stdout)
                self.assertIn("Error in %s" % s1, proc.stderr)

    def test_default_rundir(self):
        """Without -r the outputs go to the working directory, with paths built from "./"."""
        listing = self.juncfile_list(MAIN_SAMPLES)
        cwd = self.make_rundir("cwd")
        self.run_ok("-j", listing, "-m", "10", cwd=cwd)
        self.assert_outputs(main_expected("./", prefix="leafcutter"), cwd)

    def test_defaults(self):
        """-o leafcutter, -m 30, -p 0.001 and -l 100000 when not given.

        Each default sits on a boundary: the 5000:* introns have exactly 30 reads (kept)
        and the cluster at 200000 has 29 (dropped); 300000:300201 has 30/30000 = 0.001 of
        its cluster (kept) and 250:301 has 31/31031 < 0.001 (removed); 5000:105000 is
        100000 long (kept) and 5000:105001 is 100001 (dropped).
        """
        listing = self.juncfile_list([("s1.junc", [
            junc("chr1", 100, 201, 20000), junc("chr1", 100, 301, 11000),
            junc("chr1", 250, 301, 31),
            junc("chr1", 5000, 6001, 30), junc("chr1", 5000, 7001, 30),
            junc("chr1", 5000, 105000, 30), junc("chr1", 5000, 105001, 30),
            junc("chr1", 200000, 200101, 14), junc("chr1", 200000, 200201, 15),
            junc("chr1", 300000, 300101, 29970), junc("chr1", 300000, 300201, 30)])])
        self.run_ok("-j", listing, "-r", self.rundir)
        self.assertEqual(self.output("leafcutter_pooled"),
                         "chr1:+ 100:201:20000 100:301:11000 250:301:31 \n"
                         "chr1:+ 5000:6001:30 5000:7001:30 5000:105000:30 \n"
                         "chr1:+ 200000:200101:14 200000:200201:15 \n"
                         "chr1:+ 300000:300101:29970 300000:300201:30 \n")
        self.assertEqual(self.output("leafcutter_refined"),
                         "chr1:+ 100:201:20000 100:301:11000 \n"
                         "chr1:+ 5000:6001:30 5000:7001:30 5000:105000:30 \n"
                         "chr1:+ 300000:300101:29970 300000:300201:30 \n")
        self.assertEqual(self.output("leafcutter_perind.counts.gz"),
                         "chrom s1\n"
                         "chr1:100:201:clu_1_+ 20000/31000\n"
                         "chr1:100:301:clu_1_+ 11000/31000\n"
                         "chr1:5000:105000:clu_2_+ 30/90\n"
                         "chr1:5000:6001:clu_2_+ 30/90\n"
                         "chr1:5000:7001:clu_2_+ 30/90\n"
                         "chr1:300000:300101:clu_3_+ 29970/30000\n"
                         "chr1:300000:300201:clu_3_+ 30/30000\n")
        self.assertEqual(self.output("leafcutter_perind_numers.counts.gz"),
                         "s1\n"
                         "chr1:100:201:clu_1_+ 20000\n"
                         "chr1:100:301:clu_1_+ 11000\n"
                         "chr1:5000:105000:clu_2_+ 30\n"
                         "chr1:5000:6001:clu_2_+ 30\n"
                         "chr1:5000:7001:clu_2_+ 30\n"
                         "chr1:300000:300101:clu_3_+ 29970\n"
                         "chr1:300000:300201:clu_3_+ 30\n")

    def test_quiet(self):
        listing = self.juncfile_list(MAIN_SAMPLES)
        loud_dir, quiet_dir = self.make_rundir("loud"), self.make_rundir("quiet")
        loud = self.run_ok("-j", listing, "-r", loud_dir, "-o", "main", "-m", "10")
        quiet = self.run_ok("-j", listing, "-r", quiet_dir, "-o", "main", "-m", "10", "-q")
        long_dir = self.make_rundir("quiet-long")
        quiet_long = self.run_ok("-j", listing, "-r", long_dir, "-o", "main", "-m", "10",
                                 "--quiet")
        self.assertIn("scanning 1/3: %s" % os.path.join(self.indir, "s1.junc"), loud.stderr)
        self.assertIn("Sorting 1/3: ", loud.stderr)
        for proc in (quiet, quiet_long):
            self.assertNotIn("scanning", proc.stderr)
            self.assertNotIn("Sorting", proc.stderr)
        self.assert_outputs(main_expected(quiet_dir), quiet_dir)
        self.assert_outputs(main_expected(long_dir), long_dir)

    def test_maxintronlen(self):
        listing = self.juncfile_list([("s1.junc", [
            junc("chr1", 1000, 1050, 5), junc("chr1", 1000, 1100, 5),
            junc("chr1", 1000, 1101, 5)])])
        expected = {
            "100": "chrom s1\n"
                   "chr1:1000:1050:clu_1_+ 5/10\n"
                   "chr1:1000:1100:clu_1_+ 5/10\n",
            "101": "chrom s1\n"
                   "chr1:1000:1050:clu_1_+ 5/15\n"
                   "chr1:1000:1100:clu_1_+ 5/15\n"
                   "chr1:1000:1101:clu_1_+ 5/15\n",
            "50": "chrom s1\n",
        }
        for flag in ("-l", "--maxintronlen"):
            for length, counts in expected.items():
                with self.subTest(flag=flag, length=length):
                    rundir = self.make_rundir("out%s%s" % (flag, length))
                    self.run_ok("-j", listing, "-r", rundir, "-m", "3", flag, length)
                    self.assertEqual(self.output("leafcutter_perind.counts.gz", rundir),
                                     counts)
                    if length == "50":
                        self.assertEqual(
                            self.output("leafcutter_perind_numers.counts.gz", rundir), "s1\n")

    def test_minclureads(self):
        """The cluster at 100 has 20 reads: kept at -m 10, dropped at -m 30."""
        listing = self.juncfile_list([("s1.junc", [
            junc("chr1", 100, 201, 10), junc("chr1", 100, 301, 10),
            junc("chr1", 1000, 1101, 30), junc("chr1", 1000, 1201, 30)])])
        expected = {
            "10": "chrom s1\n"
                  "chr1:100:201:clu_1_+ 10/20\n"
                  "chr1:100:301:clu_1_+ 10/20\n"
                  "chr1:1000:1101:clu_2_+ 30/60\n"
                  "chr1:1000:1201:clu_2_+ 30/60\n",
            "30": "chrom s1\n"
                  "chr1:1000:1101:clu_1_+ 30/60\n"
                  "chr1:1000:1201:clu_1_+ 30/60\n",
        }
        for flag in ("-m", "--minclureads"):
            for reads, counts in expected.items():
                with self.subTest(flag=flag, reads=reads):
                    rundir = self.make_rundir("out%s%s" % (flag, reads))
                    self.run_ok("-j", listing, "-r", rundir, flag, reads)
                    self.assertEqual(self.output("leafcutter_perind.counts.gz", rundir),
                                     counts)

    def test_mincluratio(self):
        """200:500 has 4 of 124 reads (0.032): removed at -p 0.05, kept at -p 0.03."""
        listing = self.juncfile_list([("s1.junc", [
            junc("chr1", 100, 200, 30), junc("chr1", 100, 250, 30),
            junc("chr1", 200, 500, 4), junc("chr1", 400, 500, 30),
            junc("chr1", 400, 550, 30)])])
        expected = {
            "0.05": "chrom s1\n"
                    "chr1:100:200:clu_1_+ 30/60\n"
                    "chr1:100:250:clu_1_+ 30/60\n"
                    "chr1:400:500:clu_2_+ 30/60\n"
                    "chr1:400:550:clu_2_+ 30/60\n",
            "0.03": "chrom s1\n"
                    "chr1:100:200:clu_1_+ 30/124\n"
                    "chr1:100:250:clu_1_+ 30/124\n"
                    "chr1:200:500:clu_1_+ 4/124\n"
                    "chr1:400:500:clu_1_+ 30/124\n"
                    "chr1:400:550:clu_1_+ 30/124\n",
        }
        for flag in ("-p", "--mincluratio"):
            for ratio, counts in expected.items():
                with self.subTest(flag=flag, ratio=ratio):
                    rundir = self.make_rundir("out%s%s" % (flag, ratio))
                    self.run_ok("-j", listing, "-r", rundir, "-m", "3", flag, ratio)
                    self.assertEqual(self.output("leafcutter_perind.counts.gz", rundir),
                                     counts)

    def test_cluster(self):
        """Counts below 3 are used: with -c nothing is pooled or filtered."""
        listing = self.juncfile_list([
            ("s1.junc", [junc("chr1", 100, 201, 1), junc("chr1", 100, 301, 2),
                         junc("chr1", 5000, 5101, 4)]),
            ("s2.junc", [junc("chr1", 100, 201, 3), junc("chr1", 5000, 5201, 6),
                         junc("chr1", 5000, 5101, 1)])])
        clusters = os.path.join(self.indir, "clusters.txt")
        write_text(clusters, "chr1:+ 100:201:0 100:301:0 \nchr1:+ 5000:5101:0 5000:5201:0 \n")
        for flag in ("-c", "--cluster"):
            with self.subTest(flag=flag):
                rundir = self.make_rundir("out" + flag)
                self.run_ok("-j", listing, "-r", rundir, flag, clusters)
                self.assert_outputs({
                    "leafcutter_sortedlibs":
                        "%s/s1.junc.leafcutter.sorted.gz\n"
                        "%s/s2.junc.leafcutter.sorted.gz\n" % (rundir, rundir),
                    "s1.junc.leafcutter.sorted.gz":
                        "chrom s1\n"
                        "chr1:100:201:clu_1_+ 1/3\n"
                        "chr1:100:301:clu_1_+ 2/3\n"
                        "chr1:5000:5101:clu_2_+ 4/4\n"
                        "chr1:5000:5201:clu_2_+ 0/4\n",
                    "s2.junc.leafcutter.sorted.gz":
                        "chrom s2\n"
                        "chr1:100:201:clu_1_+ 3/3\n"
                        "chr1:100:301:clu_1_+ 0/3\n"
                        "chr1:5000:5101:clu_2_+ 1/7\n"
                        "chr1:5000:5201:clu_2_+ 6/7\n",
                    "leafcutter_perind.counts.gz":
                        "chrom s1 s2\n"
                        "chr1:100:201:clu_1_+ 1/3 3/3\n"
                        "chr1:100:301:clu_1_+ 2/3 0/3\n"
                        "chr1:5000:5101:clu_2_+ 4/4 1/7\n"
                        "chr1:5000:5201:clu_2_+ 0/4 6/7\n",
                    "leafcutter_perind_numers.counts.gz":
                        "s1 s2\n"
                        "chr1:100:201:clu_1_+ 1 3\n"
                        "chr1:100:301:clu_1_+ 2 0\n"
                        "chr1:5000:5101:clu_2_+ 4 1\n"
                        "chr1:5000:5201:clu_2_+ 0 6\n",
                }, rundir)

    def test_nochromcheck(self):
        listing = self.juncfile_list([("s1.junc", [
            junc("1", 100, 201, 5), junc("1", 100, 301, 5),
            junc("chrM", 100, 201, 6), junc("chrM", 100, 301, 6),
            junc("1", 100, 201, 50, "?"), junc("1", 100, 401, 50, "?")])])
        checked = "chrom s1\n1:100:201:clu_1_+ 5/10\n1:100:301:clu_1_+ 5/10\n"
        unchecked = checked + "chrM:100:201:clu_2_+ 6/12\nchrM:100:301:clu_2_+ 6/12\n"
        for label, extra, counts in (("none", [], checked),
                                     ("short", ["-k", "1"], unchecked),
                                     ("long", ["--nochromcheck", "1"], unchecked)):
            with self.subTest(label=label):
                rundir = self.make_rundir("out-" + label)
                self.run_ok("-j", listing, "-r", rundir, "-m", "3", *extra)
                self.assertEqual(self.output("leafcutter_perind.counts.gz", rundir), counts)

    def test_includeconst(self):
        """-C writes the counts under *.constcounts.gz names."""
        listing = self.juncfile_list([("s1.junc", [
            junc("chr1", 100, 201, 10), junc("chr1", 100, 301, 10),
            junc("chr1", 1000, 1101, 30), junc("chr1", 1000, 1201, 30)])])
        for flag in ("-C", "--includeconst"):
            with self.subTest(flag=flag):
                rundir = self.make_rundir("out" + flag)
                self.run_ok("-j", listing, "-r", rundir, "-m", "10", flag)
                self.assertEqual(self.output("leafcutter_perind.constcounts.gz", rundir),
                                 "chrom s1\n"
                                 "chr1:100:201:clu_1_+ 10/20\n"
                                 "chr1:100:301:clu_1_+ 10/20\n"
                                 "chr1:1000:1101:clu_2_+ 30/60\n"
                                 "chr1:1000:1201:clu_2_+ 30/60\n")
                self.assertEqual(self.output("leafcutter_perind_numers.constcounts.gz",
                                             rundir),
                                 "s1\n"
                                 "chr1:100:201:clu_1_+ 10\n"
                                 "chr1:100:301:clu_1_+ 10\n"
                                 "chr1:1000:1101:clu_2_+ 30\n"
                                 "chr1:1000:1201:clu_2_+ 30\n")
                for name in ("leafcutter_perind.counts.gz",
                             "leafcutter_perind_numers.counts.gz"):
                    self.assertFalse(os.path.exists(os.path.join(rundir, name)))

    def test_no_juncfiles_option(self):
        proc = self.run_script("-r", self.rundir)
        self.assertEqual(proc.stderr, "Error: no junction file provided...\n")
        self.assertEqual(os.listdir(self.rundir), [])

    def test_missing_junction_file(self):
        listing = self.juncfile_list([("s1.junc", [junc("chr1", 100, 201, 5)])])
        missing = os.path.join(self.indir, "missing.junc")
        with open(listing, "a") as f:
            f.write(missing + "\n")
        proc = self.run_script("-j", listing, "-r", self.rundir)
        self.assertIn("%s does not exist... check your junction files.\n" % missing,
                      proc.stderr)
        self.assertEqual(os.listdir(self.rundir), [])


if __name__ == "__main__":
    unittest.main()
