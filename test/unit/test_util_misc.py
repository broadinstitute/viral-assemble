# Unit tests for util.misc.py

__author__ = "dpark@broadinstitute.org"

import os, random, collections
import unittest
import logging
import subprocess
import util.misc
import util.file
import pytest


class TestRunAndPrint(unittest.TestCase):
    
    def testBasicRunSuccess(self):
        util.misc.run_and_print(['cat', '/dev/null'],
            silent=False, buffered=False, check=False)
        util.misc.run_and_print(['cat', '/dev/null'],
            silent=False, buffered=False, check=True)
        util.misc.run_and_print(['cat', '/dev/null'],
            silent=True, buffered=False, check=False)
        util.misc.run_and_print(['cat', '/dev/null'],
            silent=True, buffered=False, check=True)
        util.misc.run_and_print(['cat', '/dev/null'],
            silent=False, buffered=True, check=False)
        util.misc.run_and_print(['cat', '/dev/null'],
            silent=False, buffered=True, check=True)
        util.misc.run_and_print(['cat', '/dev/null'],
            silent=True, buffered=True, check=False)
        util.misc.run_and_print(['cat', '/dev/null'],
            silent=True, buffered=True, check=True)
            
    def testBasicRunFailDontCare(self):
        util.misc.run_and_print(['cat', '/notdev/notnull'],
            silent=False, buffered=False, check=False)
        util.misc.run_and_print(['cat', '/notdev/notnull'],
            silent=True, buffered=False, check=False)
        util.misc.run_and_print(['cat', '/notdev/notnull'],
            silent=False, buffered=True, check=False)
        util.misc.run_and_print(['cat', '/notdev/notnull'],
            silent=True, buffered=True, check=False)

    def testBasicRunFailAndCatch(self):
        self.assertRaises(subprocess.CalledProcessError,
            util.misc.run_and_print, ['cat', '/notdev/notnull'],
            silent=False, buffered=False, check=True)
        self.assertRaises(subprocess.CalledProcessError,
            util.misc.run_and_print, ['cat', '/notdev/notnull'],
            silent=False, buffered=True, check=True)
        self.assertRaises(subprocess.CalledProcessError,
            util.misc.run_and_print, ['cat', '/notdev/notnull'],
            silent=True, buffered=False, check=True)
        self.assertRaises(subprocess.CalledProcessError,
            util.misc.run_and_print, ['cat', '/notdev/notnull'],
            silent=True, buffered=True, check=True)


class TestFeatureSorter(unittest.TestCase):

    def testBasicSortingWithOverlap(self):
        fs = util.misc.FeatureSorter((
            ('abca', 10, 20),
            ('abca', 25, 35),
            ('abca', 15, 30),
        ))
        self.assertEqual(
            list(fs.get_features()),
            [
                ('abca', 10, 20, '+', None),
                ('abca', 15, 30, '+', None),
                ('abca', 25, 35, '+', None),
            ]
        )

    def testBasicIntervalsWithOverlap(self):
        fs = util.misc.FeatureSorter((
            ('abca', 10, 20),
            ('abca', 25, 35),
            ('abca', 15, 30),
        ))
        self.assertEqual(
            list(fs.get_intervals()),
            [
                ('abca', 10, 14, 1, [('abca', 10, 20, '+', None),]),
                ('abca', 15, 20, 2, [('abca', 10, 20, '+', None),('abca', 15, 30, '+', None),]),
                ('abca', 21, 24, 1, [('abca', 15, 30, '+', None),]),
                ('abca', 25, 30, 2, [('abca', 15, 30, '+', None),('abca', 25, 35, '+', None),]),
                ('abca', 31, 35, 1, [('abca', 25, 35, '+', None),]),
            ]
        )

    def testDisjointAndOverlappingIntervals(self):
        fs = util.misc.FeatureSorter((
            ('abca', 10, 20),
            ('abca', 80, 90, '+', None),
            ('abca', 25, 35, '-'),
            ('abca', 15, 30),
        ))
        self.assertEqual(
            list(fs.get_intervals()),
            [
                ('abca', 10, 14, 1, [('abca', 10, 20, '+', None),]),
                ('abca', 15, 20, 2, [('abca', 10, 20, '+', None),('abca', 15, 30, '+', None),]),
                ('abca', 21, 24, 1, [('abca', 15, 30, '+', None),]),
                ('abca', 25, 30, 2, [('abca', 15, 30, '+', None),('abca', 25, 35, '-', None),]),
                ('abca', 31, 35, 1, [('abca', 25, 35, '-', None),]),
                ('abca', 36, 79, 0, []),
                ('abca', 80, 90, 1, [('abca', 80, 90, '+', None),]),
            ]
        )

    def testMultiChrWindowedFeatures(self):
        fs = util.misc.FeatureSorter((
            ('abca', 10, 20),
            ('aaaa', 17, 33),
            ('abca', 80, 90),
            ('abca', 25, 35),
            ('abca', 15, 30),
        ))
        self.assertEqual(
            list(fs.get_features('abca', 11, 22)),
            [
                ('abca', 10, 20, '+', None),
                ('abca', 15, 30, '+', None),
            ]
        )

    def testOpenWindowRight(self):
        fs = util.misc.FeatureSorter((
            ('abca', 10, 20),
            ('aaaa', 17, 33),
            ('abca', 80, 90),
            ('abca', 25, 35),
            ('abca', 15, 30),
        ))
        self.assertEqual(
            list(fs.get_features('abca', 22)),
            [
                ('abca', 15, 30, '+', None),
                ('abca', 25, 35, '+', None),
                ('abca', 80, 90, '+', None),
            ]
        )

    def testOpenWindowLeft(self):
        fs = util.misc.FeatureSorter((
            ('abca', 10, 20),
            ('aaaa', 17, 33),
            ('abca', 80, 90),
            ('abca', 25, 35),
            ('abca', 15, 30),
        ))
        self.assertEqual(
            list(fs.get_features('abca', right=18)),
            [
                ('abca', 10, 20, '+', None),
                ('abca', 15, 30, '+', None),
            ]
        )

    def testMultiChrWithPayloadIntervals(self):
        fs = util.misc.FeatureSorter((
            ('abca', 10, 20),
            ('aaaa', 17, 33, '-', [100, 'name', []]),
            ('abca', 80, 90, '+', ['other info']),
            ('abca', 25, 35, '-'),
            ('abca', 15, 30),
        ))
        self.assertEqual(
            list(fs.get_intervals()),
            [
                ('abca', 10, 14, 1, [('abca', 10, 20, '+', None),]),
                ('abca', 15, 20, 2, [('abca', 10, 20, '+', None),('abca', 15, 30, '+', None),]),
                ('abca', 21, 24, 1, [('abca', 15, 30, '+', None),]),
                ('abca', 25, 30, 2, [('abca', 15, 30, '+', None),('abca', 25, 35, '-', None),]),
                ('abca', 31, 35, 1, [('abca', 25, 35, '-', None),]),
                ('abca', 36, 79, 0, []),
                ('abca', 80, 90, 1, [('abca', 80, 90, '+', ['other info']),]),
                ('aaaa', 17, 33, 1, [('aaaa', 17, 33, '-', [100, 'name', []]),]),
            ]
        )


class TestConfigIncludes(unittest.TestCase):

    def testConfigIncludes(self):

        def test_fn(f): return os.path.join(util.file.get_test_input_path(), 'TestUtilMisc', f)
        cfg1 = util.misc.load_config(test_fn('cfg1.yaml'))
        cfg2 = util.misc.load_config(test_fn('cfg2.yaml'), std_includes=[test_fn('cfg_std.yaml')],
                                     param_renamings={'std_param_A_old': 'std_param_A_new'})
        
        self.assertIn('paramA', cfg2)
        self.assertEqual(cfg2["env_vars"]["var_A"],1)
        self.assertEqual(cfg2["env_vars"]["var_B"],3)
        self.assertEqual(cfg2["env_vars"]["var_C"],4)
        with self.assertRaises(KeyError): cfg2["env_vars"]["var_E"]
        with self.assertRaises(KeyError): cfg2["var_A"]
        self.assertFalse(cfg2["paramZ"])
        self.assertTrue(cfg1["paramZ"])
        self.assertEqual(cfg2["empty_subtree"]["x"], 1)

        self.assertEqual(cfg2["std_methods"], ['a','b','c'])
        self.assertEqual(cfg2["stage1"]["stage2"]["step_num"],3)
        self.assertEqual(cfg2["stage1"]["stage2"]["step_list"],[5,10,15])
        self.assertEqual(cfg2["stage1"]["stage3"]["step_list"],[3,33])
        self.assertEqual(cfg1["stage1"]["stage3"]["step_list"],[51,101,151])

        self.assertEqual(cfg2["std_param_A_new"], 111)  # specified as std_param_A_old in cfg1.yaml

        self.assertEqual(util.misc.load_config(test_fn('empty.yaml')), {})

def test_as_type():
    """Test util.misc.as_type()"""

    as_type = util.misc.as_type

    test_data = (
        ('1', int, 1, int),
        ('1', (int,), 1, int),
        ('1', (int, float), 1, int),
        ('1', (float, int), 1., float),
        ('1.2', (float, int), 1.2, float),
        ('1.2', (int, float), 1.2, float),
        (1, int, 1, int),
        ('1.', (int, float), 1., float),
        (1., int, 1, int),
        (1.2, (int, float), 1, int),
        ('1e3', (int, float), 1000., float),
        ('1e3', (float, int), 1000., float),
        ('1.1e3', (int, float), 1100., float),
        ('-1.1e3', (float, int), -1100., float),
    )

    for val, types, out_val, out_type in test_data:
        result = as_type(val, types)
        assert result == out_val
        assert type(result) == out_type

    err_data = (
        ('1.', int), ('1.', (int,)), ('1e3', int), ('1e3', (int,)), ([1], int),
        ('', int), ('', (int, float)), ('', (float, int))
    )

    for val, types in err_data:
        with pytest.raises(TypeError):
            as_type(val, types)

@pytest.mark.parametrize("iter_d", [False, True])
@pytest.mark.parametrize("iter_subset", [False, True])
def test_subdict(iter_d, iter_subset):
    """Test util.misc.subdict()"""
    def subdict(d, subset): 
        return util.misc.subdict(iter(d.items()) if iter_d else d, 
                                 iter(subset) if iter_subset else subset)

    test_data = (
        ({}, {}, {}),
        ({1:2}, {}, {}),
        ({1:2}, {1}, {1:2}),
        ({1:2}, {2}, {}),
        ({1:2}, {1,2}, {1:2}),
        ({1:2,3:4}, {1,2,3}, {1:2,3:4}),
        ({1:2,3:4}, {2,3,5}, {3:4}),
    )

    for d, subset, expected in test_data:
        assert subdict(d, subset) == expected

        assert set(subdict(d, subset).keys()) == (set(d.keys()) & set(subset))
        assert subdict(d, {}) == {}
        assert subdict(d, d.keys()) == d
        assert subdict(d, list(d.keys())*2) == d

def test_chk():
    chk = util.misc.chk
    chk(True, 'no error')
    assert chk(1) == 1
    assert chk('1', 'nonempty string') == '1'
    with pytest.raises(RuntimeError):
        chk(False)
    with pytest.raises(RuntimeError):
        chk('', "string is empty")
    with pytest.raises(RuntimeError):
        chk(2 == 3, 'Something wrong')
    with pytest.raises(TypeError):
        chk(isinstance(None, int), 'Expected an int', TypeError)

def test_first_non_None(*args):
    first_non_None = util.misc.first_non_None
    assert first_non_None(1) == 1
    assert first_non_None(1, None) == 1
    assert first_non_None(None, False) == False
    assert first_non_None(None, None, 0) == 0
    assert first_non_None(None, '', None, True) == ''
    for v in (False, (), '', 0):
        assert first_non_None(v) == v

    with pytest.raises(ValueError):
        first_non_None()
    with pytest.raises(ValueError):
        first_non_None(None)
    with pytest.raises(ValueError):
        first_non_None(None, None)

def test_transform_json_data():
    transform_json_data = util.misc.transform_json_data
    def handle_node(val, path):
        logging.info('HANDLE_NODE: val={} path={}'.format(val, path))
        return val
    transform_json_data(1, handle_node)

def test_json_gather_leaf_jpaths():
    test_data = (
        (1, {(): 1}),
        ([1, 2], {(0,): 1, (1,): 2}),
        ({'1': 2}, {('1',): 2}),
        ({'1': 2, '3': [4, '5']}, {('1',): 2, ('3', 0): 4, ('3', 1): '5'}),
    )
    for inp, expected_out in test_data:
        assert(util.misc.json_gather_leaf_jpaths(inp) == expected_out)

def test_map_vals():
    map_vals = util.misc.map_vals
    assert tuple(map_vals({})) == ()
    assert tuple(map_vals({1:2})) == (2,)
    assert tuple(map_vals({1:2,3:4})) == (2,4)
