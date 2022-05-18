import unittest
from disgrid import data


class TestDataMethods(unittest.TestCase):
    def setUp(self):
        pass

    def test_init(self):
        self.data = data.Data("samples/example")


if __name__ == '__main__':
    unittest.main()
