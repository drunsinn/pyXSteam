import math
import re
import unittest
from collections import defaultdict

from pyXSteam import tables as table_data

TABLE_NAME = re.compile(r"^(Table(?:\d+|xx|A\d+_\d+)(?:_T3[a-z]{2})?|Table)_")


class TableDataTests(unittest.TestCase):
    def test_columns_in_each_table_have_the_same_length(self):
        for data_class in vars(table_data).values():
            if not isinstance(data_class, type) or not data_class.__module__.startswith(f"{table_data.__name__}."):
                continue

            columns_by_table = defaultdict(dict)
            for name, values in vars(data_class).items():
                match = TABLE_NAME.match(name)
                if match is None or not isinstance(values, list):
                    continue
                if any(isinstance(value, list) for value in values):
                    continue

                table_name = match.group(1)
                if data_class.__name__ == "SR2_01" and name == "Table5_n":
                    # This coefficient list is not a column of the verification data also named Table 5.
                    table_name += "_coefficients"
                columns_by_table[table_name][name] = len(values)

            for table_name, columns in columns_by_table.items():
                self.assertEqual(
                    len(set(columns.values())),
                    1,
                    f"{data_class.__name__}.{table_name} columns have inconsistent lengths: {columns}",
                )

    def test_sr4_04_coefficient_counts_match_tables_3_and_4(self):
        for table_number, expected_count in ((3, 33), (4, 35)):
            for column in ("I", "J", "n"):
                coefficients = getattr(table_data.SR4_04, f"Table{table_number}_{column}")
                self.assertEqual(
                    len(coefficients),
                    expected_count,
                    f"SR4_04.Table{table_number}_{column} has {len(coefficients)} entries; Table{table_number} specifies {expected_count}",
                )

    def test_sr5_05_coefficient_counts_match_tables_4_and_12(self):
        for table_number, counts, appendix in (
            (4, table_data.SR5_05.Table4_N, "A1"),
            (12, table_data.SR5_05.Table12_N, "A2"),
        ):
            for table_index, expected_count in enumerate(counts, start=1):
                table_name = f"Table{appendix}_{table_index}"
                for column in ("I", "J", "n"):
                    coefficients = getattr(table_data.SR5_05, f"{table_name}_{column}")
                    self.assertEqual(
                        len(coefficients),
                        expected_count,
                        f"SR5_05.{table_name}_{column} has {len(coefficients)} entries; Table{table_number} specifies {expected_count}",
                    )

    def test_sr5_05_boundary_verification_values_match(self):
        data = table_data.SR5_05
        for equation, pressure, expected in zip(data.Table3_equation, data.Table3_p, data.Table3_T):
            if equation == "T3ef":
                actual = 647.096 + 3.727888004 * (pressure - 22.064)
            else:
                suffix = equation.replace("T3", "")
                pressure = math.log(pressure) if equation in ("T3ab", "T3op") else pressure
                actual = sum(
                    coefficient * pressure**exponent
                    for coefficient, exponent in zip(
                        getattr(data, f"Table1_T3{suffix}_n"),
                        getattr(data, f"Table1_T3{suffix}_I"),
                    )
                )
            self.assertAlmostEqual(actual, expected, places=7, msg=equation)

        for equation, pressure, expected in zip(data.Table11_equation, data.Table11_p, data.Table11_T):
            pressure = math.log(pressure) if equation == "T3wx" else pressure
            actual = sum(
                coefficient * pressure**exponent
                for coefficient, exponent in zip(
                    getattr(data, f"Table9_{equation}_n"),
                    getattr(data, f"Table9_{equation}_I"),
                )
            )
            self.assertAlmostEqual(actual, expected, places=7, msg=equation)

    def test_sr5_05_volume_verification_values_match(self):
        data = table_data.SR5_05
        for prefix, params_table, samples_table in (
            ("A1", 4, 5),
            ("A2", 12, 13),
        ):
            subregions = getattr(data, f"Table{samples_table}_subregion")
            pressures = getattr(data, f"Table{samples_table}_p")
            temperatures = getattr(data, f"Table{samples_table}_T")
            expected_volumes = getattr(data, f"Table{samples_table}_v")

            for subregion, pressure, temperature, expected in zip(subregions, pressures, temperatures, expected_volumes):
                letter = subregion[-1]
                start = "a" if prefix == "A1" else "u"
                index = ord(letter) - ord(start)
                coefficients_index = index + 1
                params = [
                    getattr(data, f"Table{params_table}_{name}")[index] for name in ("v_star", "p_star", "T_star", "a", "b", "c", "d", "e")
                ]
                v_star, p_star, T_star, a, b, c, d, e = params
                exponents_i = getattr(data, f"Table{prefix}_{coefficients_index}_I")
                exponents_j = getattr(data, f"Table{prefix}_{coefficients_index}_J")
                n = getattr(data, f"Table{prefix}_{coefficients_index}_n")
                pi, theta = pressure / p_star, temperature / T_star

                if subregion == "3n":
                    omega = math.exp(
                        sum(
                            coefficient * (pi - a) ** exponent_i * (theta - b) ** exponent_j
                            for coefficient, exponent_i, exponent_j in zip(n, exponents_i, exponents_j)
                        )
                    )
                else:
                    omega = (
                        sum(
                            coefficient * (pi - a) ** (c * exponent_i) * (theta - b) ** (d * exponent_j)
                            for coefficient, exponent_i, exponent_j in zip(n, exponents_i, exponents_j)
                        )
                        ** e
                    )

                actual = v_star * omega
                self.assertTrue(
                    math.isclose(actual, expected, rel_tol=1e-8),
                    f"{subregion} at p={pressure} MPa, T={temperature} K: got {actual}, expected {expected}",
                )
