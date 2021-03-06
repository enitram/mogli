########################################################################################################################
#     mogli - molecular graph library                                                                                  #
#                                                                                                                      #
#     Copyright (C) 2016-2019  Martin S. Engler                                                                        #
#                                                                                                                      #
#     This program is free software: you can redistribute it and/or modify                                             #
#     it under the terms of the GNU Lesser General Public License as published                                         #
#     by the Free Software Foundation, either version 3 of the License, or                                             #
#     (at your option) any later version.                                                                              #
#                                                                                                                      #
#     This program is distributed in the hope that it will be useful,                                                  #
#     but WITHOUT ANY WARRANTY; without even the implied warranty of                                                   #
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the                                                     #
#     GNU General Public License for more details.                                                                     #
#                                                                                                                      #
#     You should have received a copy of the GNU Lesser General Public License                                         #
#     along with this program.  If not, see <https://www.gnu.org/licenses/>.                                           #
########################################################################################################################

import unittest


TIMEOUT = 60
TIMEOUT_BIG = 360

ETHYL = """@nodes
partial_charge label label2 atomType coordX coordY coordZ initColor
-0.048         0     C1     12       -0.765 -0.000 0.000  1
0.016          1     H3     20       -1.164 -0.813 0.619  2
0.016          2     H4     20       -1.164 -0.129 -1.013 3
0.016          3     H5     20       -1.164 0.942  0.394  4
@edges
    label
0 1 0
0 2 1
0 3 2"""

ETHANE_1 = """@nodes
partial_charge label label2 atomType    coordX coordY coordZ initColor
-0.048         0     C1     12          -0.765 -0.000 0.000  1
-0.048         1     C2     12          0.766  0.000  -0.000 2
0.016          2     H3     20          -1.164 -0.813 0.619  3
0.016          3     H4     20          -1.164 -0.129 -1.013 4
0.016          4     H5     20          -1.164 0.942  0.394  5
0.016          5     H6     20          1.164  0.129  1.013  6
0.016          6     H7     20          1.164  0.813  -0.619 7
0.016          7     H8     20          1.164  -0.942 -0.394 8
@edges
    label
0 1 0
0 2 1
0 3 2
0 4 3
1 5 4
1 6 5
1 7 6"""

ETHANE_2 = """@nodes
partial_charge label label2 atomType coordX coordY coordZ initColor
0.016          0     H3     20       -1.164 -0.813 0.619  1
0.016          1     H4     20       -1.164 -0.129 -1.013 2
0.016          2     H5     20       -1.164 0.942  0.394  3
0.016          3     H6     20       1.164  0.129  1.013  4
0.016          4     H7     20       1.164  0.813  -0.619 5
0.016          5     H8     20       1.164  -0.942 -0.394 6
-0.048         6     C1     12       -0.765 -0.000 0.000  7
-0.048         7     C2     12       0.766  0.000  -0.000 8
@edges
    label
0 6 0
1 6 1
2 6 2
3 7 3
4 7 4
5 7 5
6 7 6"""

ETHANE_PROPS = """@nodes
label label2 atomType isC    iNum  dNum
0     C1     12       1      1     1.0
1     C2     12       1      2     2.0
2     H3     20       0      3     3.0
3     H4     20       0      4     4.0
4     H5     20       0      5     5.0
5     H6     20       0      6     6.0
6     H7     20       0      7     7.0
7     H8     20       0      8     8.0
@edges
    label
0 1 0
0 2 1
0 3 2
0 4 3
1 5 4
1 6 5
1 7 6"""

GENERIC_1 = """@nodes
label atomType
0     1
1     1
2     2
3     3
4     4
5     5
6     6
7     1
8     1
@edges
    label
0 2 0
1 2 1
2 3 2
3 4 3
3 5 4
5 6 5
5 7 6
6 8 7"""

GENERIC_2 = """@nodes
label atomType
0     1
1     1
2     2
3     3
4     7
5     5
6     6
7     1
8     1
@edges
    label
0 2 0
1 2 1
2 3 2
3 4 3
3 5 4
5 6 5
5 7 6
6 8 7"""

MOLID3246 = """@nodes
partial_charge	label	label2	atomType	coordX	coordY	coordZ	initColor
-0.493	1	O3	3	-2.59	-0.833	-0.418	33
0.005	10	H1	20	-1.647	-0.04	1.25	24
0.398	11	H2	21	-1.149	-1.374	-1.763	23
0.049	12	H6	20	1.007	0.068	1.292	22
-0.725	13	OC	2	3.482	0.538	1.599	21
-0.725	14	OB	2	3.924	-2.044	1.417	20
0.928	15	P1	30	3.574	-0.732	0.773	19
-0.375	16	OA	3	2.181	-0.907	-0.109	18
0.063	17	C1	12	1.031	-0.103	0.208	17
0.061	18	H9	20	5.752	-1.06	-2.016	16
0.061	19	H8	20	5.349	-2.334	-0.83	15
-0.821	2	O33	2	-3.824	-1.379	1.803	32
0.061	20	H7	20	4.099	-1.726	-1.95	14
0.022	21	CC	12	4.974	-1.447	-1.352	13
-0.416	22	OD	3	4.652	-0.405	-0.427	12
0.356	23	C4	12	-1.5	1.356	-0.407	11
-0.046	24	H3	20	-1.551	1.269	-1.507	10
0.427	25	H4	21	-3.386	1.534	0.139	9
-0.726	26	O4	3	-2.57	2.133	0.088	8
0.023	27	H5	20	0.992	1.09	-1.581	7
-0.009	28	H4	20	-0.148	2.184	1.036	6
0.434	29	H6	21	2.678	1.649	0.534	5
-0.821	3	O32	2	-4.635	0.572	0.313	31
0.421	30	H5	21	0.782	3.613	-0.588	4
-0.736	31	O6	3	2.206	2.045	-0.245	3
0.321	32	C6	12	1.044	1.256	-0.499	2
-0.669	33	O5	3	-0.15	3.342	-0.685	1
0.152	34	C5	12	-0.191	2.057	-0.061	0
-0.821	4	O31	2	-4.883	-1.862	-0.516	30
0.934	5	P3	30	-4.128	-0.893	0.385	29
-0.631	6	O2	3	-0.212	-1.184	-1.565	28
0.057	7	C2	12	-0.234	-0.876	-0.17	27
0.171	8	C3	12	-1.523	-0.089	0.16	26
0.07	9	H2	20	-0.242	-1.814	0.408	25
@edges
		label
1	5	0
1	8	1
2	5	2
3	5	3
4	5	4
6	7	5
6	11	6
7	8	7
7	9	8
7	17	9
8	10	10
8	23	11
12	17	12
13	15	13
14	15	14
15	16	15
15	22	16
16	17	17
17	32	18
18	21	19
19	21	20
20	21	21
21	22	22
23	24	23
23	26	24
23	34	25
25	26	26
27	32	27
28	34	28
29	31	29
30	33	30
31	32	31
32	34	32
33	34	33"""

PACLITAXEL = """@nodes
partial_charge	label	label2	atomType	coordX	coordY	coordZ	initColor
0	1	X	20	0	0	0	0
0	2	X	12	0	0	0	0
0	3	X	20	0	0	0	0
0	4	X	20	0	0	0	0
0	5	X	12	0	0	0	0
0	6	X	1	0	0	0	0
0	7	X	4	0	0	0	0
0	8	X	13	0	0	0	0
0	9	X	12	0	0	0	0
0	10	X	20	0	0	0	0
0	11	X	20	0	0	0	0
0	12	X	4	0	0	0	0
0	13	X	12	0	0	0	0
0	14	X	20	0	0	0	0
0	15	X	12	0	0	0	0
0	16	X	20	0	0	0	0
0	17	X	20	0	0	0	0
0	18	X	12	0	0	0	0
0	19	X	20	0	0	0	0
0	20	X	3	0	0	0	0
0	21	X	21	0	0	0	0
0	22	X	13	0	0	0	0
0	23	X	12	0	0	0	0
0	24	X	1	0	0	0	0
0	25	X	12	0	0	0	0
0	26	X	20	0	0	0	0
0	27	X	20	0	0	0	0
0	28	X	20	0	0	0	0
0	29	X	12	0	0	0	0
0	30	X	20	0	0	0	0
0	31	X	12	0	0	0	0
0	32	X	20	0	0	0	0
0	33	X	4	0	0	0	0
0	34	X	12	0	0	0	0
0	35	X	1	0	0	0	0
0	36	X	12	0	0	0	0
0	37	X	12	0	0	0	0
0	38	X	20	0	0	0	0
0	39	X	12	0	0	0	0
0	40	X	20	0	0	0	0
0	41	X	12	0	0	0	0
0	42	X	20	0	0	0	0
0	43	X	12	0	0	0	0
0	44	X	20	0	0	0	0
0	45	X	12	0	0	0	0
0	46	X	20	0	0	0	0
0	47	X	13	0	0	0	0
0	48	X	3	0	0	0	0
0	49	X	21	0	0	0	0
0	50	X	12	0	0	0	0
0	51	X	20	0	0	0	0
0	52	X	20	0	0	0	0
0	53	X	13	0	0	0	0
0	54	X	12	0	0	0	0
0	55	X	12	0	0	0	0
0	56	X	20	0	0	0	0
0	57	X	20	0	0	0	0
0	58	X	20	0	0	0	0
0	59	X	12	0	0	0	0
0	60	X	12	0	0	0	0
0	61	X	20	0	0	0	0
0	62	X	20	0	0	0	0
0	63	X	20	0	0	0	0
0	64	X	12	0	0	0	0
0	65	X	20	0	0	0	0
0	66	X	4	0	0	0	0
0	67	X	12	0	0	0	0
0	68	X	1	0	0	0	0
0	69	X	12	0	0	0	0
0	70	X	20	0	0	0	0
0	71	X	3	0	0	0	0
0	72	X	21	0	0	0	0
0	73	X	12	0	0	0	0
0	74	X	20	0	0	0	0
0	75	X	12	0	0	0	0
0	76	X	12	0	0	0	0
0	77	X	20	0	0	0	0
0	78	X	6	0	0	0	0
0	79	X	21	0	0	0	0
0	80	X	12	0	0	0	0
0	81	X	1	0	0	0	0
0	82	X	12	0	0	0	0
0	83	X	12	0	0	0	0
0	84	X	20	0	0	0	0
0	85	X	12	0	0	0	0
0	86	X	20	0	0	0	0
0	87	X	12	0	0	0	0
0	88	X	20	0	0	0	0
0	89	X	12	0	0	0	0
0	90	X	20	0	0	0	0
0	91	X	12	0	0	0	0
0	92	X	20	0	0	0	0
0	93	X	12	0	0	0	0
0	94	X	20	0	0	0	0
0	95	X	12	0	0	0	0
0	96	X	20	0	0	0	0
0	97	X	12	0	0	0	0
0	98	X	20	0	0	0	0
0	99	X	12	0	0	0	0
0	100	X	20	0	0	0	0
0	101	X	12	0	0	0	0
0	102	X	20	0	0	0	0
0	103	X	4	0	0	0	0
0	104	X	12	0	0	0	0
0	105	X	1	0	0	0	0
0	106	X	12	0	0	0	0
0	107	X	20	0	0	0	0
0	108	X	20	0	0	0	0
0	109	X	20	0	0	0	0
0	110	X	12	0	0	0	0
0	111	X	20	0	0	0	0
0	112	X	20	0	0	0	0
0	113	X	20	0	0	0	0
@edges
		label
1	2	0
2	3	1
2	4	2
2	5	3
5	6	4
5	7	5
7	8	6
8	9	7
8	13	8
8	29	9
9	10	10
9	11	11
9	12	12
12	13	13
13	14	14
13	15	15
15	16	16
15	17	17
15	18	18
18	19	19
18	20	20
18	22	21
20	21	22
22	23	23
22	25	24
22	29	25
23	24	26
23	101	27
25	26	28
25	27	29
25	28	30
29	30	31
29	31	32
31	32	33
31	33	34
31	47	35
33	34	36
34	35	37
34	36	38
36	37	39
36	45	40
37	38	41
37	39	42
39	40	43
39	41	44
41	42	45
41	43	46
43	44	47
43	45	48
45	46	49
47	48	50
47	50	51
47	53	52
48	49	53
50	51	54
50	52	55
50	64	56
53	54	57
53	55	58
53	110	59
54	59	60
54	101	61
55	56	62
55	57	63
55	58	64
59	60	65
59	64	66
60	61	67
60	62	68
60	63	69
64	65	70
64	66	71
66	67	72
67	68	73
67	69	74
69	70	75
69	71	76
69	73	77
71	72	78
73	74	79
73	75	80
73	78	81
75	76	82
75	93	83
76	77	84
76	99	85
78	79	86
78	80	87
80	81	88
80	82	89
82	83	90
82	87	91
83	84	92
83	85	93
85	86	94
85	91	95
87	88	96
87	89	97
89	90	98
89	91	99
91	92	100
93	94	101
93	95	102
95	96	103
95	97	104
97	98	105
97	99	106
99	100	107
101	102	108
101	103	109
103	104	110
104	105	111
104	106	112
106	107	113
106	108	114
106	109	115
110	111	116
110	112	117
110	113	118"""


class TestIO(unittest.TestCase):

    def test_read_lgf(self):
        from mogli import Molecule

        mol = Molecule()
        mol.read_lgf(ETHANE_1)

        self.assertEqual(mol.get_atom_count(), 8)

        c1 = mol.get_node_by_id(0)
        c2 = mol.get_node_by_id(1)

        props = mol.get_properties()

        self.assertIn('label2', props)
        self.assertEqual(mol.get_property(c1, 'label2'), 'C1')
        self.assertEqual(mol.get_property(c2, 'label2'), 'C2')

        c1_edges = {(0, 1), (0, 2), (0, 3), (0, 4)}
        c2_edges = {(0, 1), (1, 5), (1, 6), (1, 7)}

        edges = c1_edges | c2_edges

        for e in mol.get_inc_edge_iter(c1):
            id1 = mol.get_id(mol.get_u(e))
            id2 = mol.get_id(mol.get_v(e))
            edge = (id1, id2) if id1 < id2 else (id2, id1)
            self.assertIn(edge, c1_edges)

        for e in mol.get_inc_edge_iter(c2):
            id1 = mol.get_id(mol.get_u(e))
            id2 = mol.get_id(mol.get_v(e))
            edge = (id1, id2) if id1 < id2 else (id2, id1)
            self.assertIn(edge, c2_edges)

        for e in mol.get_edge_iter():
            id1 = mol.get_id(mol.get_u(e))
            id2 = mol.get_id(mol.get_v(e))
            edge = (id1, id2) if id1 < id2 else (id2, id1)
            self.assertIn(edge, edges)

    def test_read_lgf_properties(self):
        from mogli import LGFIOConfig, Molecule

        config = LGFIOConfig("label", "atomType")
        config.add_string_node_prop("label2")\
            .add_bool_node_prop("isC")\
            .add_int_node_prop("iNum")\
            .add_double_node_prop("dNum")

        mol = Molecule()
        mol.read_lgf(ETHANE_PROPS, config)

        self.assertEqual(mol.get_atom_count(), 8)

        props = mol.get_properties()
        self.assertEqual(len(props), 4)

        self.assertIn("label2", props)
        self.assertIn("isC", props)
        self.assertIn("iNum", props)
        self.assertIn("dNum", props)

        v = mol.get_node_by_id(0)
        self.assertEqual(mol.get_property(v, "label2"), "C1")
        self.assertTrue(mol.get_property(v, "isC"))
        self.assertEqual(mol.get_property(v, "iNum"), 1)
        self.assertEqual(mol.get_property(v, "dNum"), 1.0)

        v = mol.get_node_by_id(1)
        self.assertEqual(mol.get_property(v, "label2"), "C2")
        self.assertTrue(mol.get_property(v, "isC"))
        self.assertEqual(mol.get_property(v, "iNum"), 2)
        self.assertEqual(mol.get_property(v, "dNum"), 2.0)

        v = mol.get_node_by_id(2)
        self.assertEqual(mol.get_property(v, "label2"), "H3")
        self.assertFalse(mol.get_property(v, "isC"))
        self.assertEqual(mol.get_property(v, "iNum"), 3)
        self.assertEqual(mol.get_property(v, "dNum"), 3.0)

        v = mol.get_node_by_id(3)
        self.assertEqual(mol.get_property(v, "label2"), "H4")
        self.assertFalse(mol.get_property(v, "isC"))
        self.assertEqual(mol.get_property(v, "iNum"), 4)
        self.assertEqual(mol.get_property(v, "dNum"), 4.0)

        v = mol.get_node_by_id(4)
        self.assertEqual(mol.get_property(v, "label2"), "H5")
        self.assertFalse(mol.get_property(v, "isC"))
        self.assertEqual(mol.get_property(v, "iNum"), 5)
        self.assertEqual(mol.get_property(v, "dNum"), 5.0)

        v = mol.get_node_by_id(5)
        self.assertEqual(mol.get_property(v, "label2"), "H6")
        self.assertFalse(mol.get_property(v, "isC"))
        self.assertEqual(mol.get_property(v, "iNum"), 6)
        self.assertEqual(mol.get_property(v, "dNum"), 6.0)

        v = mol.get_node_by_id(6)
        self.assertEqual(mol.get_property(v, "label2"), "H7")
        self.assertFalse(mol.get_property(v, "isC"))
        self.assertEqual(mol.get_property(v, "iNum"), 7)
        self.assertEqual(mol.get_property(v, "dNum"), 7.0)

        v = mol.get_node_by_id(7)
        self.assertEqual(mol.get_property(v, "label2"), "H8")
        self.assertFalse(mol.get_property(v, "isC"))
        self.assertEqual(mol.get_property(v, "iNum"), 8)
        self.assertEqual(mol.get_property(v, "dNum"), 8.0)


class TestIsomorphism(unittest.TestCase):

    def test_canonization(self):
        from mogli import Molecule, LGFIOConfig, Canonization

        mol1 = Molecule()
        mol2 = Molecule()
        mol3 = Molecule()
        
        config = LGFIOConfig("label", "atomType")
        mol1.read_lgf(ETHANE_1)
        mol2.read_lgf(ETHANE_2)
        mol3.read_lgf(GENERIC_1, config)

        c1 = Canonization(mol1)
        c2 = Canonization(mol2)
        c3 = Canonization(mol3)

        self.assertEqual(c1.get_canonization(), c2.get_canonization())
        self.assertEqual(c1.get_colors(), c2.get_colors())
        self.assertNotEqual(c1.get_canonization(), c3.get_canonization())
        self.assertNotEqual(c1.get_colors(), c3.get_colors())


    def test_fcanonization(self):
        from mogli import Molecule, LGFIOConfig, GenerationType, maximal_common_fragments, FragmentCanonization

        mol1 = Molecule()
        mol2 = Molecule()
        mol3 = Molecule()
        mol4 = Molecule()

        config = LGFIOConfig("label", "atomType")
        mol1.read_lgf(ETHANE_1)
        mol2.read_lgf(ETHANE_1)
        mol3.read_lgf(ETHANE_1)
        mol4.read_lgf(ETHYL, config)

        t1, frag1, _, _ = maximal_common_fragments(mol1, mol2, 1, TIMEOUT, GenerationType.NO_OPT)
        t2, frag2, _ ,_ = maximal_common_fragments(mol1, mol3, 1, TIMEOUT, GenerationType.NO_OPT)
        t3, frag3, _ ,_ = maximal_common_fragments(mol1, mol4, 1, TIMEOUT, GenerationType.NO_OPT)

        self.assertTrue(t1)
        self.assertTrue(t2)
        self.assertTrue(t3)

        self.assertTrue(len(frag1) > 0)
        self.assertTrue(len(frag2) > 0)
        self.assertTrue(len(frag3) > 0)

        frag1 = sorted(frag1, key=lambda x: x.get_atom_count(), reverse=True)
        frag2 = sorted(frag2, key=lambda x: x.get_atom_count(), reverse=True)
        frag3 = sorted(frag3, key=lambda x: x.get_atom_count(), reverse=True)

        f1 = FragmentCanonization(frag1[0])
        f2 = FragmentCanonization(frag2[0])
        f3 = FragmentCanonization(frag3[0])

        self.assertEqual(f1.get_canonization(), f2.get_canonization())
        self.assertEqual(f1.get_colors(), f2.get_colors())
        self.assertNotEqual(f1.get_canonization(), f3.get_canonization())
        self.assertNotEqual(f1.get_colors(), f3.get_colors())

    def test_isomorphism(self):
        from mogli import Molecule

        mol1 = Molecule()
        mol2 = Molecule()
        mol3 = Molecule()

        mol1.read_lgf(ETHANE_1)
        mol2.read_lgf(ETHANE_2)
        mol3.read_lgf(ETHYL)

        self.assertTrue(mol1.is_isomorphic(mol2))
        self.assertFalse(mol1.is_isomorphic(mol3))

    def test_subgraph_isomorphism(self):
        from mogli import Molecule, are_subgraph_isomorphic

        mol1 = Molecule()
        mol2 = Molecule()

        mol1.read_lgf(ETHYL)
        mol2.read_lgf(ETHANE_2)

        self.assertTrue(are_subgraph_isomorphic(mol1, mol2)[0])
        self.assertFalse(are_subgraph_isomorphic(mol2, mol1)[0])


class TestMatching(unittest.TestCase):

    def test_mcf_isomorphic_graphs(self):
        from mogli import Molecule, maximal_common_fragments, GenerationType

        mol1, mol2 = Molecule(), Molecule()

        mol1.read_lgf(ETHANE_1)
        mol2.read_lgf(ETHANE_2)

        t1, frag_noopt, _, _ = maximal_common_fragments(mol1, mol2, 1, TIMEOUT, GenerationType.NO_OPT)
        t2, frag_deg1, _, _ = maximal_common_fragments(mol1, mol2, 1, TIMEOUT, GenerationType.DEG_1)
        t3, frag_uncon, _, _ = maximal_common_fragments(mol1, mol2, 1, TIMEOUT, GenerationType.UNCON)
        t4, frag_uncon_deg1, _, _ = maximal_common_fragments(mol1, mol2, 1, TIMEOUT, GenerationType.UNCON_DEG_1)

        self.assertTrue(t1)
        self.assertTrue(t2)
        self.assertTrue(t3)
        self.assertTrue(t4)

        self.assertTrue(len(frag_noopt) > 0)
        self.assertTrue(len(frag_deg1) > 0)
        self.assertTrue(len(frag_uncon) > 0)
        self.assertTrue(len(frag_uncon_deg1) > 0)

        frag_noopt = sorted(frag_noopt, key=lambda x: x.get_atom_count(), reverse=True)
        frag_deg1 = sorted(frag_deg1, key=lambda x: x.get_atom_count(), reverse=True)
        frag_uncon = sorted(frag_uncon, key=lambda x: x.get_atom_count(), reverse=True)
        frag_uncon_deg1 = sorted(frag_uncon_deg1, key=lambda x: x.get_atom_count(), reverse=True)

        self.assertEqual(frag_noopt[0].get_atom_count(), mol1.get_atom_count())
        self.assertEqual(frag_deg1[0].get_atom_count(), mol1.get_atom_count())
        self.assertEqual(frag_uncon[0].get_atom_count(), mol1.get_atom_count())
        self.assertEqual(frag_uncon_deg1[0].get_atom_count(), mol1.get_atom_count())

    def test_mcf_isomorphic_graphs_max(self):
        from mogli import Molecule, maximal_common_fragments, GenerationType

        mol1, mol2 = Molecule(), Molecule()

        mol1.read_lgf(ETHANE_1)
        mol2.read_lgf(ETHANE_2)

        t1, frag_noopt, _, _ = maximal_common_fragments(mol1, mol2, 1, TIMEOUT, GenerationType.NO_OPT, maximum = True)
        t2, frag_deg1, _, _ = maximal_common_fragments(mol1, mol2, 1, TIMEOUT, GenerationType.DEG_1, maximum =  True)
        t3, frag_uncon, _, _ = maximal_common_fragments(mol1, mol2, 1, TIMEOUT, GenerationType.UNCON, maximum =  True)
        t4, frag_uncon_deg1, _, _ = maximal_common_fragments(mol1, mol2, 1, TIMEOUT, GenerationType.UNCON_DEG_1, maximum = True)

        self.assertTrue(t1)
        self.assertTrue(t2)
        self.assertTrue(t3)
        self.assertTrue(t4)

        self.assertTrue(len(frag_noopt) > 0)
        self.assertTrue(len(frag_deg1) > 0)
        self.assertTrue(len(frag_uncon) > 0)
        self.assertTrue(len(frag_uncon_deg1) > 0)

        for f in frag_noopt:
            self.assertEqual(f.get_atom_count(), mol1.get_atom_count())
        for f in frag_deg1:
            self.assertEqual(f.get_atom_count(), mol1.get_atom_count())
        for f in frag_uncon:
            self.assertEqual(f.get_atom_count(), mol1.get_atom_count())
        for f in frag_uncon_deg1:
            self.assertEqual(f.get_atom_count(), mol1.get_atom_count())

    def test_mcf_isomorphic_graphs_big(self):
        from mogli import Molecule, maximal_common_fragments, GenerationType

        mol1, mol2 = Molecule(), Molecule()

        mol1.read_lgf(MOLID3246)
        mol2.read_lgf(MOLID3246)

        t1, frag_noopt, _, _ = maximal_common_fragments(mol1, mol2, 1, TIMEOUT_BIG, GenerationType.NO_OPT)
        t2, frag_deg1, _, _ = maximal_common_fragments(mol1, mol2, 1, TIMEOUT_BIG, GenerationType.DEG_1)
        t3, frag_uncon_deg1, _, _ = maximal_common_fragments(mol1, mol2, 1, TIMEOUT_BIG, GenerationType.UNCON_DEG_1)

        self.assertTrue(t1)
        self.assertTrue(t2)
        self.assertTrue(t3)

        self.assertTrue(len(frag_noopt) > 0)
        self.assertTrue(len(frag_deg1) > 0)
        self.assertTrue(len(frag_uncon_deg1) > 0)

        frag_noopt = sorted(frag_noopt, key=lambda x: x.get_atom_count(), reverse=True)
        frag_deg1 = sorted(frag_deg1, key=lambda x: x.get_atom_count(), reverse=True)
        frag_uncon_deg1 = sorted(frag_uncon_deg1, key=lambda x: x.get_atom_count(), reverse=True)

        self.assertEqual(frag_noopt[0].get_atom_count(), mol1.get_atom_count())
        self.assertEqual(frag_deg1[0].get_atom_count(), mol1.get_atom_count())
        self.assertEqual(frag_uncon_deg1[0].get_atom_count(), mol1.get_atom_count())

    def test_mcf_isomorphic_graphs_large(self):
        from mogli import Molecule, maximal_common_fragments, GenerationType

        mol1, mol2 = Molecule(), Molecule()

        mol1.read_lgf(PACLITAXEL)
        mol2.read_lgf(PACLITAXEL)

        t1, _, _, _ = maximal_common_fragments(mol1, mol2, 1, TIMEOUT, GenerationType.NO_OPT)
        t2, frag_uncon_deg1, _, _ = maximal_common_fragments(mol1, mol2, 1, TIMEOUT_BIG, GenerationType.UNCON_DEG_1, maximum = True)

        self.assertFalse(t1)
        self.assertTrue(t2)

        self.assertTrue(len(frag_uncon_deg1) > 0)

        frag_uncon_deg1 = sorted(frag_uncon_deg1, key=lambda x: x.get_atom_count(), reverse=True)

        self.assertEqual(frag_uncon_deg1[0].get_atom_count(), mol1.get_atom_count())

    def test_mcf_subisomorphic_graphs_1(self):
        from mogli import Molecule, maximal_common_fragments, GenerationType

        mol1, mol2 = Molecule(), Molecule()

        mol1.read_lgf(ETHYL)
        mol2.read_lgf(ETHANE_1)

        t1, frag_noopt, _, _ = maximal_common_fragments(mol1, mol2, 0, TIMEOUT, GenerationType.NO_OPT)
        t2, frag_deg1, _, _ = maximal_common_fragments(mol1, mol2, 0, TIMEOUT, GenerationType.DEG_1)
        t3, frag_uncon, _, _ = maximal_common_fragments(mol1, mol2, 0, TIMEOUT, GenerationType.UNCON)
        t4, frag_uncon_deg1, _, _ = maximal_common_fragments(mol1, mol2, 0, TIMEOUT, GenerationType.UNCON_DEG_1)

        self.assertTrue(t1)
        self.assertTrue(t2)
        self.assertTrue(t3)
        self.assertTrue(t4)

        self.assertTrue(len(frag_noopt) > 0)
        self.assertTrue(len(frag_deg1) > 0)
        self.assertTrue(len(frag_uncon) > 0)
        self.assertTrue(len(frag_uncon_deg1) > 0)

        frag_noopt = sorted(frag_noopt, key=lambda x: x.get_atom_count(), reverse=True)
        frag_deg1 = sorted(frag_deg1, key=lambda x: x.get_atom_count(), reverse=True)
        frag_uncon = sorted(frag_uncon, key=lambda x: x.get_atom_count(), reverse=True)
        frag_uncon_deg1 = sorted(frag_uncon_deg1, key=lambda x: x.get_atom_count(), reverse=True)

        self.assertEqual(frag_noopt[0].get_atom_count(), mol1.get_atom_count())
        self.assertEqual(frag_deg1[0].get_atom_count(), mol1.get_atom_count())
        self.assertEqual(frag_uncon[0].get_atom_count(), mol1.get_atom_count())
        self.assertEqual(frag_uncon_deg1[0].get_atom_count(), mol1.get_atom_count())

    def test_mcf_subisomorphic_graphs_2(self):
        from mogli import Molecule, maximal_common_fragments, GenerationType

        mol1, mol2 = Molecule(), Molecule()

        mol1.read_lgf(ETHYL)
        mol2.read_lgf(ETHANE_2)

        t1, frag_noopt, _, _ = maximal_common_fragments(mol1, mol2, 1, TIMEOUT, GenerationType.NO_OPT)
        t2, frag_deg1, _, _ = maximal_common_fragments(mol1, mol2, 1, TIMEOUT, GenerationType.DEG_1)
        t3, frag_uncon, _, _ = maximal_common_fragments(mol1, mol2, 1, TIMEOUT, GenerationType.UNCON)
        t4, frag_uncon_deg1, _, _ = maximal_common_fragments(mol1, mol2, 1, TIMEOUT, GenerationType.UNCON_DEG_1)

        self.assertTrue(t1)
        self.assertTrue(t2)
        self.assertTrue(t3)
        self.assertTrue(t4)

        self.assertTrue(len(frag_noopt) > 0)
        self.assertTrue(len(frag_deg1) > 0)
        self.assertTrue(len(frag_uncon) > 0)
        self.assertTrue(len(frag_uncon_deg1) > 0)

        frag_noopt = sorted(frag_noopt, key=lambda x: x.get_atom_count(), reverse=True)
        frag_deg1 = sorted(frag_deg1, key=lambda x: x.get_atom_count(), reverse=True)
        frag_uncon = sorted(frag_uncon, key=lambda x: x.get_atom_count(), reverse=True)
        frag_uncon_deg1 = sorted(frag_uncon_deg1, key=lambda x: x.get_atom_count(), reverse=True)

        self.assertEqual(frag_noopt[0].get_atom_count(), 2)
        self.assertEqual(frag_deg1[0].get_atom_count(), 2)
        self.assertEqual(frag_uncon[0].get_atom_count(), 2)
        self.assertEqual(frag_uncon_deg1[0].get_atom_count(), 2)

    def test_mcf_graphs_no_match(self):
        from mogli import Molecule, LGFIOConfig, maximal_common_fragments, GenerationType

        mol1, mol2 = Molecule(), Molecule()

        config = LGFIOConfig('label', 'atomType')
        mol1.read_lgf(ETHANE_1)
        mol2.read_lgf(GENERIC_1, config)

        t1, frag_noopt, _, _ = maximal_common_fragments(mol1, mol2, 1, TIMEOUT, GenerationType.NO_OPT)
        t2, frag_deg1, _, _ = maximal_common_fragments(mol1, mol2, 1, TIMEOUT, GenerationType.DEG_1)
        t3, frag_uncon, _, _ = maximal_common_fragments(mol1, mol2, 1, TIMEOUT, GenerationType.UNCON)
        t4, frag_uncon_deg1, _, _ = maximal_common_fragments(mol1, mol2, 1, TIMEOUT, GenerationType.UNCON_DEG_1)

        self.assertTrue(t1)
        self.assertTrue(t2)
        self.assertTrue(t3)
        self.assertTrue(t4)

        self.assertEqual(len(frag_noopt), 0)
        self.assertEqual(len(frag_deg1), 0)
        self.assertEqual(len(frag_uncon), 0)
        self.assertEqual(len(frag_uncon_deg1), 0)

    def test_mcf_custom_periodic_table(self):
        import copy
        from mogli import Molecule, PeriodicTable, LGFIOConfig, maximal_common_fragments, GenerationType

        custom_table = PeriodicTable(PeriodicTable.get_default()).make_equivalent(4, 7)
        mol1, mol2, mol3, mol4 = Molecule(), Molecule(), Molecule(custom_table), Molecule(custom_table)

        config = LGFIOConfig('label', 'atomType')
        mol1.read_lgf(GENERIC_1, config)
        mol2.read_lgf(GENERIC_2, config)
        mol3.read_lgf(GENERIC_1, config)
        mol4.read_lgf(GENERIC_2, config)

        t1, frag_default, _, _ = maximal_common_fragments(mol1, mol2, 1, TIMEOUT, GenerationType.UNCON_DEG_1, True)
        t2, frag_custom1, _, _ = maximal_common_fragments(mol3, mol4, 1, TIMEOUT, GenerationType.NO_OPT, True)
        t3, frag_custom2, _, _ = maximal_common_fragments(mol3, mol4, 1, TIMEOUT, GenerationType.DEG_1, True)
        t4, frag_custom3, _, _ = maximal_common_fragments(mol3, mol4, 1, TIMEOUT, GenerationType.UNCON, True)
        t5, frag_custom4, _, _ = maximal_common_fragments(mol3, mol4, 1, TIMEOUT, GenerationType.UNCON_DEG_1, True)

        self.assertTrue(t1)
        self.assertTrue(t2)
        self.assertTrue(t3)
        self.assertTrue(t4)
        self.assertTrue(t5)

        self.assertTrue(len(frag_default) > 0)
        self.assertTrue(len(frag_custom1) > 0)
        self.assertTrue(len(frag_custom2) > 0)
        self.assertTrue(len(frag_custom3) > 0)
        self.assertTrue(len(frag_custom4) > 0)

        self.assertEqual(frag_default[0].get_atom_count(), 5)
        self.assertEqual(frag_custom1[0].get_atom_count(), mol1.get_atom_count())
        self.assertEqual(frag_custom2[0].get_atom_count(), mol1.get_atom_count())
        self.assertEqual(frag_custom3[0].get_atom_count(), mol1.get_atom_count())
        self.assertEqual(frag_custom4[0].get_atom_count(), mol1.get_atom_count())

    def test_atomic_fragments(self):
        from mogli import Molecule, atomic_fragments

        mol = Molecule()
        mol.read_lgf(ETHANE_1)

        frag0, _ = atomic_fragments(mol, 0)
        frag1, _ = atomic_fragments(mol, 1)
        frag2, _ = atomic_fragments(mol, 2)

        self.assertEqual(len(frag0), 8)
        self.assertEqual(len(frag1), 8)
        self.assertEqual(len(frag2), 8)


class TestPacking(unittest.TestCase):

    def test_hash_canonization(self):
        from mogli import Molecule, LGFIOConfig, Canonization, hash_canonization

        mol1 = Molecule()
        mol2 = Molecule()
        mol3 = Molecule()

        config = LGFIOConfig("label", "atomType")
        mol1.read_lgf(ETHANE_1)
        mol2.read_lgf(ETHANE_2)
        mol3.read_lgf(GENERIC_1, config)

        c1 = hash_canonization(Canonization(mol1))
        c2 = hash_canonization(Canonization(mol2))
        c3 = hash_canonization(Canonization(mol3))

        self.assertNotEqual(c1, c2)
        self.assertNotEqual(c1, c3)

    def test_hash_fcanonization(self):
        from mogli import Molecule, LGFIOConfig, GenerationType, maximal_common_fragments,\
            FragmentCanonization, hash_fcanonization

        mol1 = Molecule()
        mol2 = Molecule()
        mol3 = Molecule()
        mol4 = Molecule()

        config = LGFIOConfig("label", "atomType")
        mol1.read_lgf(ETHANE_1)
        mol2.read_lgf(ETHANE_1)
        mol3.read_lgf(ETHANE_1)
        mol4.read_lgf(ETHYL, config)

        t1, frag1, _, _ = maximal_common_fragments(mol1, mol2, 1, TIMEOUT, GenerationType.NO_OPT)
        t2, frag2, _ ,_ = maximal_common_fragments(mol1, mol3, 1, TIMEOUT, GenerationType.NO_OPT)
        t3, frag3, _ ,_ = maximal_common_fragments(mol1, mol4, 1, TIMEOUT, GenerationType.NO_OPT)

        self.assertTrue(t1)
        self.assertTrue(t2)
        self.assertTrue(t3)

        self.assertTrue(len(frag1) > 0)
        self.assertTrue(len(frag2) > 0)
        self.assertTrue(len(frag3) > 0)

        frag1 = sorted(frag1, key=lambda x: x.get_atom_count(), reverse=True)
        frag2 = sorted(frag2, key=lambda x: x.get_atom_count(), reverse=True)
        frag3 = sorted(frag3, key=lambda x: x.get_atom_count(), reverse=True)

        f1 = hash_fcanonization(FragmentCanonization(frag1[0]))
        f2 = hash_fcanonization(FragmentCanonization(frag2[0]))
        f3 = hash_fcanonization(FragmentCanonization(frag3[0]))

        self.assertEqual(f1, f2)
        self.assertNotEqual(f1, f3)

    def test_pack_canonization(self):
        from mogli import Molecule, Canonization, pack_canonization, unpack_canonization

        mol = Molecule()
        mol.read_lgf(ETHANE_1)

        c1 = Canonization(mol)
        c2 = unpack_canonization(pack_canonization(c1))

        self.assertEqual(c1.get_colors(), c2.get_colors())
        self.assertEqual(c1.get_canonization(), c2.get_canonization())
        self.assertEqual(c1.get_node_order(), c2.get_node_order())


    def test_pack_fcanonization(self):
        from mogli import Molecule, GenerationType, maximal_common_fragments, \
            FragmentCanonization, pack_fcanonization, unpack_fcanonization

        mol1 = Molecule()
        mol2 = Molecule()

        mol1.read_lgf(ETHANE_1)
        mol2.read_lgf(ETHANE_2)

        t1, frag, _, _ = maximal_common_fragments(mol1, mol2, 1, TIMEOUT, GenerationType.NO_OPT)

        self.assertTrue(t1)

        frag = sorted(frag, key=lambda x: x.get_atom_count(), reverse=True)

        f1 = FragmentCanonization(frag[0])
        foo = pack_fcanonization(f1)
        f2 = unpack_fcanonization(foo)

        self.assertEqual(f1.get_colors(), f2.get_colors())
        self.assertEqual(f1.get_canonization(), f2.get_canonization())
        self.assertEqual(f1.get_node_order(), f2.get_node_order())
        self.assertEqual(f1.get_core_nodes(), f2.get_core_nodes())

    def test_pack_fragment(self):
        from mogli import Molecule, GenerationType, maximal_common_fragments, pack_fragment, unpack_fragment

        mol1 = Molecule()
        mol2 = Molecule()

        mol1.read_lgf(ETHANE_1)
        mol2.read_lgf(ETHANE_2)

        t1, frag, _, _ = maximal_common_fragments(mol1, mol2, 1, TIMEOUT, GenerationType.NO_OPT)

        self.assertTrue(t1)

        frag = sorted(frag, key=lambda x: x.get_atom_count(), reverse=True)

        f1 = frag[0]
        f2 = unpack_fragment(pack_fragment(f1))

        self.assertEqual(f1.get_atom_count(), f2.get_atom_count())
        self.assertEqual(f1.get_core_atom_count(), f2.get_core_atom_count())

        def s(e):
            id1 = f1.get_id(f1.get_u(e))
            id2 = f1.get_id(f1.get_v(e))
            return (id1, id2) if id1 < id2 else (id2, id1)

        edges1 = sorted([s(e) for e in f1.get_edge_iter()])
        edges2 = sorted([s(e) for e in f2.get_edge_iter()])

        self.assertEqual(edges1, edges2)

        for v1 in f1.get_node_iter():
            v2 = f2.get_node_by_id(f1.get_id(v1))
            self.assertEqual(f1.is_core(v1), f2.is_core(v2))

    def test_pack_match(self):
        from mogli import Match, pack_match, unpack_match

        m1 = Match()

        m1.add_frag_to_mol(0, 1)
        m1.add_frag_to_mol(1, 2)

        m2 = unpack_match(pack_match(m1))

        i1 = m1.get_atom_ids()
        i2 = m2.get_atom_ids()

        self.assertEqual(i1, i2)
        self.assertEqual(m1.frag_to_mol(0), m2.frag_to_mol(0))
        self.assertEqual(m1.frag_to_mol(1), m2.frag_to_mol(1))
        self.assertEqual(m1.frag_to_mol(2), m2.frag_to_mol(2))

    def test_pack_molecule(self):
        from mogli import Molecule, LGFIOConfig, pack_molecule, unpack_molecule

        mol1 = Molecule()
        mol2 = Molecule()

        config = LGFIOConfig('label', 'atomType')
        config.add_string_node_prop("label2")\
            .add_bool_node_prop("isC")\
            .add_int_node_prop("iNum")\
            .add_double_node_prop("dNum")


        mol1.read_lgf(ETHANE_PROPS, config)
        mol2 = unpack_molecule(pack_molecule(mol1))

        self.assertEqual(mol1.get_atom_count(), mol2.get_atom_count())

        def s(e):
            id1 = mol1.get_id(mol1.get_u(e))
            id2 = mol1.get_id(mol1.get_v(e))
            return (id1, id2) if id1 < id2 else (id2, id1)

        edges1 = sorted([s(e) for e in mol1.get_edge_iter()])
        edges2 = sorted([s(e) for e in mol2.get_edge_iter()])

        self.assertEqual(edges1, edges2)

        props = mol2.get_properties()

        self.assertIn('label2', props)
        self.assertIn('isC', props)
        self.assertIn('iNum', props)
        self.assertIn('dNum', props)

        for v1 in mol1.get_node_iter():
            v2 = mol2.get_node_by_id(mol1.get_id(v1))
            self.assertEqual(mol1.get_property(v1, 'label2'), mol2.get_property(v2, 'label2'))
            self.assertEqual(mol1.get_property(v1, 'isC'), mol2.get_property(v2, 'isC'))
            self.assertEqual(mol1.get_property(v1, 'iNum'), mol2.get_property(v2, 'iNum'))
            self.assertEqual(mol1.get_property(v1, 'dNum'), mol2.get_property(v2, 'dNum'))


if __name__ == '__main__':
    unittest.main()
