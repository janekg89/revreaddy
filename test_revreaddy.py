import unittest
import revreaddy as rdy

class TestRevreaddyWrapper(unittest.TestCase):
	def test_properties(self):
		sim = rdy.Sim() # default implementation

		sim.timestep = 0.4
		self.assertEqual(sim.timestep, 0.4)
		with self.assertRaises(Exception):
			sim.timestep = 0.

		sim.kt = 35.
		self.assertEqual(sim.kt, 35.)
		with self.assertRaises(Exception):
			sim.kt = -1.

		sim.is_periodic = False
		self.assertEqual(sim.is_periodic, False)
		sim.is_periodic = True
		self.assertEqual(sim.is_periodic, True)

		sim.boxsize = 2.
		self.assertEqual(sim.boxsize, 2.)
		with self.assertRaises(Exception):
			sim.boxsize = -1.

		sim.use_neighborlist = False
		self.assertEqual(sim.use_neighborlist, False)
		sim.use_neighborlist = True
		self.assertEqual(sim.use_neighborlist, True)


if __name__ == '__main__':
	unittest.main()