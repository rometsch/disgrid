import unittest
from simdata import alias

class TestAliasMethods(unittest.TestCase):
    def setUp(self):
        self.a = alias.Alias()

    def test_register(self):
        self.a.register("foo", "bar")
        self.a.register("foo", {"type a" : "bar", "type b" : "baz"})

    def test_register_fail(self):
        with self.assertRaises(ValueError):
            self.a.register("foo", 2)
        with self.assertRaises(ValueError):
            self.a.register("foo", { "a" : 1, "b" : "bar"})

    def test_register_dict(self):
        self.a.register_dict({"foo" : "bar", "batz" : "boo"})
        self.assertEqual( self.a("batz"), "boo" )

    def test_match(self):
        self.a.register("foo", "bar")
        self.a.register("zoo", {"type a" : "animalpark", "type b" : "forest"})
        self.assertEqual( self.a("foo"), "bar" )
        self.assertEqual( self.a("zoo", variable_type="type a"), "animalpark")


if __name__ == '__main__':
    unittest.main()
