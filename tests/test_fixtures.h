////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//    mogli - molecular graph library                                                                                 //
//                                                                                                                    //
//    Copyright (C) 2016-2019  Martin S. Engler                                                                       //
//                                                                                                                    //
//    This program is free software: you can redistribute it and/or modify                                            //
//    it under the terms of the GNU Lesser General Public License as published                                        //
//    by the Free Software Foundation, either version 3 of the License, or                                            //
//    (at your option) any later version.                                                                             //
//                                                                                                                    //
//    This program is distributed in the hope that it will be useful,                                                 //
//    but WITHOUT ANY WARRANTY; without even the implied warranty of                                                  //
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the                                                    //
//    GNU General Public License for more details.                                                                    //
//                                                                                                                    //
//    You should have received a copy of the GNU Lesser General Public License                                        //
//    along with this program.  If not, see <https://www.gnu.org/licenses/>.                                          //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef MOGLI_TEST_FIXTURES_H
#define MOGLI_TEST_FIXTURES_H

#include <string>

static const std::string ETHYL = R"(@nodes
partial_charge	label	label2	atomType	coordX	coordY	coordZ	initColor
-0.048          0     C1      12        -0.765  -0.000  0.000   1
0.016           1     H3      20        -1.164  -0.813  0.619   2
0.016           2     H4      20        -1.164  -0.129  -1.013  3
0.016           3     H5      20        -1.164  0.942   0.394   4
@edges
		label
0 1 0
0	2	1
0	3	2)";

static const std::string ETHANE_1 = R"(@nodes
partial_charge	label	label2	atomType	coordX	coordY	coordZ	initColor
-0.048          0     C1      12        -0.765  -0.000  0.000   1
-0.048          1     C2      12        0.766   0.000   -0.000  2
0.016           2     H3      20        -1.164  -0.813  0.619   3
0.016           3     H4      20        -1.164  -0.129  -1.013  4
0.016           4     H5      20        -1.164  0.942   0.394   5
0.016           5     H6      20        1.164   0.129   1.013   6
0.016           6     H7      20        1.164   0.813   -0.619  7
0.016           7     H8      20        1.164   -0.942  -0.394  8
@edges
		label
0 1 0
0	2	1
0	3	2
0	4	3
1	5	4
1	6	5
1	7	6)";

static const std::string ETHANE_2 = R"(@nodes
partial_charge	label	label2	atomType	coordX	coordY	coordZ	initColor
0.016           0     H3      20        -1.164  -0.813  0.619   1
0.016           1     H4      20        -1.164  -0.129  -1.013  2
0.016           2     H5      20        -1.164  0.942   0.394   3
0.016           3     H6      20        1.164   0.129   1.013   4
0.016           4     H7      20        1.164   0.813   -0.619  5
0.016           5     H8      20        1.164   -0.942  -0.394  6
-0.048          6     C1      12        -0.765  -0.000  0.000   7
-0.048          7     C2      12        0.766   0.000   -0.000  8
@edges
		label
0 6 0
1 6 1
2 6 2
3 7 3
4 7 4
5 7 5
6 7 6)";

static const std::string ETHANE_PROPS = R"(@nodes
label	label2	atomType	isC	   iNum  dNum
0     C1      12        1      1     1.0
1     C2      12        1      2     2.0
2     H3      20        0      3     3.0
3     H4      20        0      4     4.0
4     H5      20        0      5     5.0
5     H6      20        0      6     6.0
6     H7      20        0      7     7.0
7     H8      20        0      8     8.0
@edges
		label
0 1 0
0	2	1
0	3	2
0	4	3
1	5	4
1	6	5
1	7	6)";

static const std::string GENERIC_1 = R"(@nodes
label	atomType
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
6 8 7
)";

static const std::string GENERIC_2 = R"(@nodes
label	atomType
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
6 8 7
)";

static const std::string PACLITAXEL = R"(@nodes
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
110	113	118)";

#endif //MOGLI_TEST_FIXTURES_H