"""
Tests the functionality in the chemistry.unit package.
"""
from __future__ import division

from chemistry import unit as u
import math
import unittest
from utils import has_numpy, numpy as np

class TestUnits(unittest.TestCase):

    def setUp(self):
        self.furlong = u.BaseUnit(u.length_dimension, "furlong", "fur")

    def assertAlmostEqualQuantities(self, item1, item2):
        try:
            val1 = item1.value_in_unit(item1.unit)
            val2 = item2.value_in_unit(item1.unit)
        except TypeError:
            raise self.failureException('Incompatible units %s and %s' %
                                        (item1.unit, item2.unit))
        self.assertAlmostEqual(val1, val2)

    def testBaseUnit(self):
        """ Tests the creation of a base unit furlong """
        furlong = u.BaseUnit(u.length_dimension, "furlong", "fur")
        furlong.define_conversion_factor_to(u.meter_base_unit, 201.168)
        self.assertEqual(furlong.conversion_factor_to(u.angstrom_base_unit),
                         201.168e10)

    def testIsUnit(self):
        """ Tests the unit.is_unit introspective function """
        self.assertTrue(u.is_unit(u.meters))
        self.assertFalse(u.is_unit(None))
        self.assertFalse(u.is_unit(10))
        # Distinguish between "units" and "quantities"
        self.assertFalse(u.is_unit(1*u.meters))

    def testCalorieConversion(self):
        """ Tests unit conversion with calories to joules """
        c = 1.0 * u.calories
        j = 1.0 * u.joules
        j2c = j.in_units_of(u.calories)
        self.assertNotEqual(c, 1.0) # Units do not equal scalars!
        self.assertIs(c.unit, u.calories)
        self.assertEqual(c._value, 1)
        self.assertEqual(u.calorie.conversion_factor_to(u.joule), 4.184)
        self.assertEqual(u.joule.conversion_factor_to(u.calorie), 1/4.184)
        self.assertIs(j2c.unit, u.calories)
        self.assertEqual(j2c._value, 1/4.184)
        self.assertEqual(c / u.calories, 1.0) # Dimensionless now
        self.assertEqual(j / u.calories, j2c._value)
        c2 = c**2
        cc = c * c
        self.assertEqual(cc, c2)
        self.assertTrue(c2.unit, u.calories**2)
        self.assertEqual(c2.in_units_of(u.joules**2), (4.184*u.joules)**2)

    def testScaledUnit(self):
        """ Tests ScaledUnit class with kilokelvins """
        kK = u.ScaledUnit(1000.0, u.kelvin, "kilokelvin", "kK")
        self.assertIs(kK.master, u.kelvin)
        self.assertEqual(kK.factor, 1000)
        self.assertEqual(str(kK), 'kilokelvin')
        self.assertEqual(kK.name, 'kilokelvin')
        self.assertEqual(kK.symbol, 'kK')

    def testComparisonOperators(self):
        """ Tests unit comparison operators """
        self.assertGreater(u.meters, u.centimeters)
        self.assertLess(u.angstroms, u.centimeters)

    def testIllegalComparison(self):
        """ Tests illegal comparisons between incompatible units """
        self.assertRaises(TypeError, lambda: u.meters > u.liters)

    def testUnitDivision(self):
        """ Tests the division of units and is_dimensionless """
        mps = u.meter / u.second
        mpm = u.meter / u.meter
        mpc = u.meter / u.centimeter
        self.assertTrue(u.is_unit(mps))
        self.assertTrue(u.is_unit(mpm))
        self.assertTrue(u.is_unit(mpc))
        self.assertFalse(mps.is_dimensionless())
        self.assertTrue(mpm.is_dimensionless())
        self.assertTrue(mpc.is_dimensionless())

    def testCompositeUnits(self):
        """ Tests the creation of a composite unit """
        mps = u.Unit({u.meter_base_unit : 1.0, u.second_base_unit : -1.0})
        self.assertTrue(u.is_unit(mps))
        self.assertEqual(str(mps), 'meter/second')

    def testUnitSystem(self):
        """ Tests the creation of a UnitSystem and its behavior """
        us = u.UnitSystem([u.ScaledUnit(1.0, u.coulomb/u.second, 'ampere', 'A'),
                           u.second_base_unit])
        self.assertEqual(us.express_unit(u.second), u.second)
        self.assertEqual(us.express_unit(u.hour), u.second)
        self.assertNotEqual(u.hour, u.second)
        self.assertEqual(us.express_unit(u.coulomb/u.second), u.ampere)
        self.assertEqual(us.express_unit(u.meter/u.second), u.meter/u.second)
        self.assertEqual(us.express_unit(u.kilometer/u.hour),
                         u.kilometer/u.second)

    def testUnitSqrt(self):
        """ Tests taking the square root of units """
        self.assertEqual((u.meter * u.meter).sqrt(), u.meter)
        self.assertEqual((u.meter**4).sqrt(), u.meter**2)

    def testUnitBadSqrt(self):
        """ Tests that illegal sqrt calls on incompatible units fails """
        mps2 = u.meters/u.second**2
        self.assertRaises(ArithmeticError, lambda: u.meter.sqrt())
        self.assertRaises(ArithmeticError, lambda: (u.meters**3).sqrt())
        self.assertRaises(ArithmeticError, lambda: mps2.sqrt())

    def testBaseScaleMix(self):
        """ Test mixing of BaseUnit and ScaledUnit instances """
        kgj = u.kilogram * u.joule
        self.assertTrue(u.is_unit(kgj))
        self.assertEqual(str(kgj.sqrt()), 'kilogram*meter/second')

    def testGetUnitAttributes(self):
        """ Tests the unit attribute `getters' """
        self.assertEqual(u.newton.get_name(), 'newton')
        self.assertEqual(u.newtons.get_name(), 'newton')
        self.assertEqual(u.newton.get_symbol(), 'N')
        self.assertEqual(u.ampere.get_symbol(), 'A')
        self.assertEqual(u.meter.get_name(), 'meter')
        self.assertEqual(u.meter.get_symbol(), 'm')

    def testPresetUnitSystems(self):
        """ Tests some of the pre-set UnitSystem's """
        self.assertEqual(u.angstrom.in_unit_system(u.si_unit_system), u.meter)
        self.assertEqual(u.angstrom.in_unit_system(u.cgs_unit_system),
                         u.centimeter)
        self.assertEqual(u.angstrom.in_unit_system(u.md_unit_system),
                         u.nanometer)
        mps = u.meter / u.second**2
        self.assertEqual(str(mps), 'meter/(second**2)')
        self.assertEqual(mps.in_unit_system(u.si_unit_system), mps)
        self.assertEqual(mps.in_unit_system(u.cgs_unit_system),
                         u.centimeter/u.second**2)
        self.assertEqual(mps.in_unit_system(u.md_unit_system),
                         u.nanometer/u.picosecond**2)

    def testIsCompatible(self):
        """ Tests the is_compatible attribute of units """
        self.assertTrue(u.meter.is_compatible(u.centimeter))
        self.assertTrue(u.centimeter.is_compatible(u.meter))
        self.assertTrue(u.meter.is_compatible(u.meter))
        self.assertTrue(u.joule.is_compatible(u.calorie))
        self.assertFalse(u.meter.is_compatible(u.kelvin))
        self.assertFalse(u.meter.is_compatible(u.meter/u.second))
        self.assertFalse(u.meter.is_compatible(u.joule))

    def testConversionFactorTo(self):
        """ Tests the "conversion_factor_to" attribute of Units """
        self.assertEqual(u.meter.conversion_factor_to(u.centimeter), 100)
        self.assertEqual(u.kilocalorie.conversion_factor_to(u.joule), 4184)
        self.assertEqual(u.calorie.conversion_factor_to(u.kilojoule), 4.184e-3)
        self.assertEqual((u.kilocalorie/u.mole/u.angstrom).conversion_factor_to(
                                u.kilojoule/u.mole/u.nanometer), 41.84)

    def testUnitString(self):
        """ Test the Unit->str casting functionality """
        # Always alphabetical order
        self.assertEqual(str(u.meter * u.second * u.second * u.kilogram),
                         'kilogram*meter*second**2')
        self.assertEqual(str(u.meter / u.second / u.second / u.kilogram),
                         'meter/(kilogram*second**2)')
        self.assertEqual(str(u.meter**3), 'meter**3')

    def testToBaseUnit(self):
        """ Tests the "get_conversion_factor_to_base_units" method """
        self.assertEqual(u.meter.get_conversion_factor_to_base_units(), 1)
        self.assertEqual(u.calorie.get_conversion_factor_to_base_units(), 4.184)
        kcpma = u.md_kilocalorie/u.mole/u.angstrom
        self.assertEqual(kcpma.get_conversion_factor_to_base_units(), 4.184)

    def testScalarQuantityMultiplyDivide(self):
        """ Tests creating a scalar Quantity object by * or / by a Unit """
        self.assertTrue(u.is_quantity(5 * u.centimeters))
        self.assertTrue(u.is_quantity(1 / u.centimeters))
        self.assertTrue(u.is_quantity(10 * u.centimeters))
        self.assertTrue(u.is_quantity(9.81 * u.meters / u.second**2))

    def testScalarQuantityConstructor(self):
        """ Tests creating a Quantity using the Quantity constructor """
        self.assertTrue(u.is_quantity(u.Quantity(5, u.centimeters)))
        self.assertTrue(u.is_quantity(u.Quantity(5, u.centimeters**-1)))

    def testValueInUnit(self):
        """ Tests the value_in_unit functionality for Quantity """
        i = 5 * u.centimeters
        self.assertEqual(i.value_in_unit(u.millimeters), 50)
        self.assertEqual(i / u.millimeters, 50)

    def testCollectionQuantities(self):
        """ Tests the use of collections as Quantity values """
        s = [1, 2, 3] * u.centimeters
        self.assertEqual(str(s), '[1, 2, 3] cm')
        self.assertTrue(u.is_quantity(s))
        s2 = s / u.millimeters
        self.assertEqual(s2, [10.0, 20.0, 30.0])
        self.assertEqual(s2, s.value_in_unit(u.millimeters))
        # Test 2-D list
        s = [[1, 2, 3], [4, 5, 6]]
        s *= u.centimeters
        self.assertTrue(u.is_quantity(s))
        s2 = s / u.millimeters
        self.assertEqual(s2, [[10.0, 20.0, 30.0], [40.0, 50.0, 60.0]])
        self.assertEqual(s.value_in_unit(u.millimeters), s2)
        # Test tuples
        s = (1, 2, 3) * u.centimeters
        self.assertTrue(u.is_quantity(s))
        self.assertEqual(str(s), '(1, 2, 3) cm')
        s2 = s / u.millimeters
        self.assertEqual(s2, (10, 20, 30))
        self.assertIsInstance(s2, tuple)
        self.assertEqual(s.value_in_unit(u.millimeters), s2)
        self.assertIsInstance(s.value_in_unit(u.millimeters), tuple)

    def testCollectionQuantityOperations(self):
        """ Tests that Quantity collections behave correctly """
        # Tests that __getitem__ returns a unit
        s = [1, 2, 3, 4] * u.angstroms
        self.assertTrue(u.is_quantity(s[0]))
        for i, val in enumerate(s):
            self.assertTrue(u.is_quantity(val))
            self.assertEqual(val, (i+1) * u.angstroms)
        # Tests that __setitem__ fails when an incompatible type is added
        def fail(s): s[0] = 5
        self.assertRaises(AttributeError, lambda: fail(s))
        def fail(s): s[0] = 5 * u.joules
        self.assertRaises(TypeError, lambda: fail(s))
        def fail(s): s[0] /= 10 * u.meters
        self.assertRaises(AttributeError, lambda: fail(s))
        # Tests that __setitem__ converts to the unit of the container
        s[0] = 1 * u.nanometers
        self.assertEqual(s[0]._value, 10)

    def testMutableQuantityOperations(self):
        " Tests that mutable Quantity objects do not get unexpectedly changed "
        # This used to be a bug -- t and s._value were the same object, so
        # changing t would also change s silently
        s = [1, 2, 3, 4] * u.angstroms
        t = s / u.angstroms
        self.assertEqual(t, [1, 2, 3, 4])
        self.assertEqual(s, u.Quantity([1, 2, 3, 4], u.angstroms))
        t[0] = 2
        self.assertEqual(t, [2, 2, 3, 4])
        self.assertEqual(s, u.Quantity([1, 2, 3, 4], u.angstroms))
        t = s.value_in_unit(u.angstroms)
        self.assertEqual(t, [1, 2, 3, 4])
        self.assertEqual(s, u.Quantity([1, 2, 3, 4], u.angstroms))
        t[0] = 2
        self.assertEqual(t, [2, 2, 3, 4])
        self.assertEqual(s, u.Quantity([1, 2, 3, 4], u.angstroms))

    def testQuantityMaths(self):
        """ Tests dimensional analysis & maths on and b/w Quantity objects """
        x = 1.3 * u.meters
        y = 75.2 * u.centimeters
        self.assertEqual((x + y) / u.meters, 2.052)
        self.assertEqual((x - y) / u.meters, 0.548)
        self.assertEqual(x / y, 1.3 / 0.752)
        self.assertEqual(x * y, 1.3*0.752*u.meters**2)
        # Now check that the left operand determines the unit
        d1 = 2.0*u.meters
        d2 = 2.0*u.nanometers
        self.assertEqual(d1 + d2, (2+2e-9)*u.meters)
        self.assertAlmostEqual((d2+d1-(2e9+2)*u.nanometers)._value, 0, places=6)
        self.assertEqual(d1 + d1, 4.0*u.meters)
        self.assertEqual(d1 - d1, 0.0*u.meters)
        self.assertEqual(d1 / d1, 1.0)
        self.assertEqual(d1 * u.meters, 2.0*u.meters**2)
        self.assertEqual(u.kilograms*(d1/u.seconds)*(d1/u.seconds),
                         4*u.kilograms*u.meters**2/u.seconds**2)
        self.assertEqual(u.kilograms*(d1/u.seconds)**2,
                         4*u.kilograms*u.meters**2/u.seconds**2)
        self.assertEqual(d1**3, 8.0*u.meters**3)
        x = d1**(3/2)
        self.assertAlmostEqual(x._value, math.sqrt(2)**3)
        self.assertEqual(x.unit, u.meters**(3/2))
        self.assertAlmostEqual((d1**0.5)._value, math.sqrt(2))
        self.assertEqual((d1**0.5).unit, u.meters**0.5)
        comp = (3.0 + 4.0j) * u.meters
        self.assertTrue(u.is_quantity(comp))
        self.assertEqual(comp.unit, u.meters)
        self.assertEqual(str(comp), '(3+4j) m')
        self.assertEqual(comp + comp, (6.0 + 8.0j)*u.meters)
        self.assertEqual(comp - comp, 0*u.meters)
        self.assertEqual(comp * comp, (3.0 + 4.0j)**2 * u.meters**2)
        self.assertAlmostEqual(comp / comp, 1)

    def testQuantityComplicatedMaths(self):
        """ Tests a complicated set of mathematical operations on a Quantity """
        s1 = 2.0
        x1 = 2
        x2 = 4.0 / 3.0
        u1 = u.kilogram * u.meter / u.second**2
        u2 = u1 * u.meter
        q1 = 1.0 * u1
        q2 = 2.0 * u2
        self.assertEqual(s1, 2.0)
        self.assertEqual(x1, 2)
        self.assertAlmostEqual(x2, 1.33333333333333)
        self.assertEqual(str(u1), 'kilogram*meter/(second**2)')
        self.assertEqual(str(u2), 'kilogram*meter**2/(second**2)')
        self.assertEqual(str(q1), '1.0 kg m/(s**2)')
        self.assertEqual(str(q2), '2.0 kg m**2/(s**2)')
        self.assertEqual(str(u1*s1), '2.0 kg m/(s**2)')
        self.assertEqual(str(s1*u1), '2.0 kg m/(s**2)')
        self.assertEqual(str(u1/s1), '0.5 kg m/(s**2)')
        self.assertEqual(str(s1/u1), '2.0 s**2/(kg m)')
        self.assertEqual(str(u1*u1), 'kilogram**2*meter**2/(second**4)')
        self.assertEqual(u1/u1, u.dimensionless)
        self.assertEqual(str(u1/u1), 'dimensionless')
        self.assertEqual(str(u1*u2), 'kilogram**2*meter**3/(second**4)')
        self.assertEqual(u1/u2, u.meters**-1)
        self.assertEqual(str(u1/u2), '/meter')
        self.assertEqual(u1**x1, u.kilogram**2*u.meter**2/(u.second**4))
        self.assertEqual(u1**(1/x1), u.kilogram**0.5*u.meter**0.5/u.second)
        self.assertEqual(u1**x2,
                         u.kilogram**(1+1/3)*u.meter**(1+1/3)/u.second**(2+2/3))

    def testQuantityComparisons(self):
        """ Tests binary comparison operators between Quantity """
        l1 = 1.0 * u.meters
        l2 = 2.0 * u.meters
        l3 = 1.0 * u.meters
        self.assertEqual(l1, l3)
        self.assertNotEqual(l1, l2)
        self.assertLessEqual(l1, l2)
        self.assertLess(l1, l2)
        self.assertGreater(l2, l1)
        self.assertGreaterEqual(l2, l1)
        self.assertLessEqual(l1, l3)
        self.assertGreaterEqual(l1, l3)

    def testChemistryProblems(self):
        """ Tests some gen-chem applications with Quantity's """
        def work(f, dx):
            return f * dx

        F = 1.0 * u.kilogram * u.meter / u.second**2
        dx = 1.0 * u.meter
        self.assertEqual(work(F, dx), 1.0 * u.joule)
        self.assertEqual(F, 1.0 * u.newton)

        def ideal_gas_law(P, V, T):
            R = u.MOLAR_GAS_CONSTANT_R
            return (P * V / (R * T)).in_units_of(u.mole)

        T = (273.0 + 37.0) * u.kelvin
        P = (1.01325e5) * u.pascals
        r = 0.5e-6 * u.meters
        V = 4/3 * math.pi * r**3
        n = ideal_gas_law(P, V, T)
        val = 4/3*math.pi*0.5e-6**3*1
        self.assertAlmostEqualQuantities(P*V, val * u.atmospheres*u.meters**3)
        self.assertAlmostEqualQuantities(n, 2.05834818672e-17 * u.mole)
        self.assertAlmostEqualQuantities(V, 5.2359833333333e-19 * u.meters**3)
        self.assertEqual(str(T), '310.0 K')
        self.assertEqual(str(u.MOLAR_GAS_CONSTANT_R), '8.31447247122 J/(K mol)')
        self.assertTrue(u.is_quantity(V))

    def testAngleQuantities(self):
        """ Tests angle measurements """
        self.assertEqual(1.0*u.radians / u.degrees, 180 / math.pi)
        self.assertTrue(u.is_quantity(1.0*u.radians))
        self.assertTrue(u.is_quantity(1.0*u.degrees))
        self.assertEqual((1.0*u.radians).in_units_of(u.degrees),
                         (180 / math.pi)*u.degrees)

    def testBadQuantityMaths(self):
        """ Tests that Maths of incompatible units properly fails """
        self.assertRaises(TypeError, lambda:1*u.meters + 1*u.liters)
        self.assertRaises(AttributeError, lambda: 1*u.liters + 5)

class TestNumpyUnits(unittest.TestCase):

    def testNumpyQuantity(self):
        """ Tests that numpy arrays can form Quantity values """
        q = u.Quantity(np.array([1, 2, 3]), u.centimeters)
        self.assertTrue(u.is_quantity(q))
        self.assertIsInstance(q._value, np.ndarray)
        self.assertTrue(np.all(q / u.millimeters == np.array([1, 2, 3])*10))
