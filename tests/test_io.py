from tonic.io import config_type, isfloat, isint, isscalar


def test_config_type_int():
    assert config_type('1') == 1


def test_config_type_float():
    assert config_type('1.75') == 1.75


def test_config_type_bool():
    assert config_type('True')


def test_config_type():
    x = [True, None, 'abc', 1, 1.5]
    assert config_type(", ".join(map(str, x))) == x


def test_isfloat():
    assert isfloat(4.3)
    assert isfloat('4.3')
    assert isfloat(4)
    assert isfloat('4')
    assert not isfloat('four')


def test_isint():
    assert isint(4)
    assert isint('4')
    assert not isint(4.3)
    assert not isint('4.3')


def test_isscalar():
    assert not isscalar([0, 1])
    assert not isscalar(('a', 'b'))
    assert isscalar(1)
    assert not isscalar('str')
