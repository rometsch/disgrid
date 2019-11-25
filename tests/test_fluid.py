import unittest
from simdata import fluid


class TestFluidMethods(unittest.TestCase):
    def setUp(self):
        self.fluid = fluid.Fluid("test")

    def test_register_wrong_type(self):
        with self.assertRaises(ValueError):
            self.fluid.register_variable("foo", "bar_type", lambda: "baz")

    def test_register_correct_type(self):
        self.fluid.register_variable("foo", "scalar", lambda: "baz")
        self.fluid.register_variable("foo", "2d", lambda: "baz")

    def test_register_non_callable(self):
        with self.assertRaises(TypeError):
            self.fluid.register_variable("foo", "2d", "baz")

    def test_register_callable(self):
        self.fluid.register_variable("foo", "2d", lambda: "baz")

    # def test_pass_wrong_varialbe_loaders_dict(self):
    #     with self.assertRaises(ValueError):
    #         fluid.Fluid("test", {"field": {}, "WRONG": {}})


if __name__ == '__main__':
    unittest.main()
