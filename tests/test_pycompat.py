from tonic.pycompat import pyrange, pyzip, iteritems


def test_pyrange():
    x = pyrange(5)
    assert list(x) == [0, 1, 2, 3, 4]


def test_pyzip():
    a = [1, 2, 3]
    b = ['a', 'b', 'c']

    z = pyzip(a, b)
    assert hasattr(z, '__iter__')
    lz = list(pyzip(a, b))
    assert lz[0] == (1, 'a')
    assert len(lz) == len(a)  # == len(b)


def test_iteritems():
    d = {'a': 0, 'b': 1, 'c': 2}
    for k, v in iteritems(d):
        assert type(k) == str
        assert type(v) == int
