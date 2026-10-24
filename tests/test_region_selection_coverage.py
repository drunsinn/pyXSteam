import unittest
from unittest.mock import patch

from pyXSteam.Constants import DiagramRegion
from pyXSteam.RegionSelection import (
    select_region_hs,
    select_region_ph,
    select_region_prho,
    select_region_ps,
    select_region_pT,
)
from pyXSteam.Regions import Region1, Region2, Region4, Region5


class RegionSelectionCoverageTester(unittest.TestCase):
    """Check region dispatch using IF97 verification states."""

    def test_select_region_pT(self):
        cases = (
            (3.0, 300.0, DiagramRegion.R1),
            (0.0035, 300.0, DiagramRegion.R2),
            (25.583702, 650.0, DiagramRegion.R3),
            (1.0, Region4.T4_p(1.0), DiagramRegion.R4),
            (0.5, 1500.0, DiagramRegion.R5),
            (0.0001, 300.0, DiagramRegion.NILL),
        )
        for pressure, temperature, expected_region in cases:
            with self.subTest(pressure=pressure, temperature=temperature):
                self.assertEqual(select_region_pT(pressure, temperature), expected_region)

    def test_select_region_ph(self):
        cases = (
            (3.0, 500.0, DiagramRegion.R1),
            (0.0035, 3000.0, DiagramRegion.R2),
            (20.0, 1700.0, DiagramRegion.R3),
            (1.0, 1500.0, DiagramRegion.R4),
            (0.5, 5219.76855, DiagramRegion.R5),
            (100.0, 6000.0, DiagramRegion.NILL),
            (0.0001, 300.0, DiagramRegion.NILL),
        )
        for pressure, enthalpy, expected_region in cases:
            with self.subTest(pressure=pressure, enthalpy=enthalpy):
                self.assertEqual(select_region_ph(pressure, enthalpy), expected_region)

    def test_select_region_ps(self):
        cases = (
            (3.0, 0.5, DiagramRegion.R1),
            (0.1, 7.5, DiagramRegion.R2),
            (20.0, 3.8, DiagramRegion.R3),
            (1.0, 5.0, DiagramRegion.R4),
            (0.5, 9.65408875, DiagramRegion.R5),
            (101.0, 5.0, DiagramRegion.NILL),
            (0.1, -1.0, DiagramRegion.NILL),
        )
        for pressure, entropy, expected_region in cases:
            with self.subTest(pressure=pressure, entropy=entropy):
                self.assertEqual(select_region_ps(pressure, entropy), expected_region)

    def test_select_region_ps_reuses_saturation_temperature(self):
        with patch(
            "pyXSteam.RegionSelection.Region4.T4_p",
            wraps=Region4.T4_p,
        ) as temperature_lookup:
            self.assertEqual(select_region_ps(1.0, 5.0), DiagramRegion.R4)

        temperature_lookup.assert_called_once_with(1.0)

    def test_select_region_hs(self):
        cases = (
            (115.331273, 0.392294792, DiagramRegion.R1),
            (2800.0, 6.5, DiagramRegion.R2),
            (1863.43019, 4.05427273, DiagramRegion.R3),
            (1800.0, 5.3, DiagramRegion.R4),
            (5219.76855, 9.65408875, DiagramRegion.NILL),
            (500.0, -1.0, DiagramRegion.NILL),
        )
        for enthalpy, entropy, expected_region in cases:
            with self.subTest(enthalpy=enthalpy, entropy=entropy):
                self.assertEqual(select_region_hs(enthalpy, entropy), expected_region)

    def test_select_region_prho(self):
        saturation_temperature = Region4.T4_p(1.0)
        liquid_volume = Region1.v1_pT(1.0, saturation_temperature)
        vapour_volume = Region2.v2_pT(1.0, saturation_temperature)
        mixture_density = 1.0 / (0.5 * liquid_volume + 0.5 * vapour_volume)

        cases = (
            (3.0, 997.8529401, DiagramRegion.R1),
            (0.0035, 0.0253219774, DiagramRegion.R2),
            (25.583702, 500.0, DiagramRegion.R3),
            (1.0, mixture_density, DiagramRegion.R4),
            (0.5, 1.0 / Region5.v5_pT(0.5, 1500.0), DiagramRegion.R5),
            (101.0, 500.0, DiagramRegion.NILL),
        )
        for pressure, density, expected_region in cases:
            with self.subTest(pressure=pressure, density=density):
                self.assertEqual(select_region_prho(pressure, density), expected_region)

    def test_select_region_prho_reuses_saturation_temperature(self):
        saturation_temperature = Region4.T4_p(1.0)
        liquid_volume = Region1.v1_pT(1.0, saturation_temperature)
        vapour_volume = Region2.v2_pT(1.0, saturation_temperature)
        mixture_density = 1.0 / (0.5 * liquid_volume + 0.5 * vapour_volume)

        with patch(
            "pyXSteam.RegionSelection.Region4.T4_p",
            wraps=Region4.T4_p,
        ) as temperature_lookup:
            self.assertEqual(select_region_prho(1.0, mixture_density), DiagramRegion.R4)

        temperature_lookup.assert_called_once_with(1.0)
