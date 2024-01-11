#ifndef eopdut_h
#define eopdut_h

/*
 * The values for UT1-UTC since 20190101
 * Created with eopdut.py on 2019-07-20T09:46:03.992188
 * Note: the DUT precision drops dramatically in around 20 days after 20190130!
 *       Please update the data in time.
 * utc (mjd), ut1-utc (sec),  nr, year,  mon,  day
*/
//const int papi::eopdutver = 20190130;
//const int papi::neopdut = 397;
double papi::eopdut[397][2] = {
  {58484.00, -0.0361533}, //    0 2019  1  1
  {58485.00, -0.0370347}, //    1 2019  1  2
  {58486.00, -0.0377546}, //    2 2019  1  3
  {58487.00, -0.0382730}, //    3 2019  1  4
  {58488.00, -0.0387241}, //    4 2019  1  5
  {58489.00, -0.0391168}, //    5 2019  1  6
  {58490.00, -0.0394338}, //    6 2019  1  7
  {58491.00, -0.0397247}, //    7 2019  1  8
  {58492.00, -0.0400892}, //    8 2019  1  9
  {58493.00, -0.0405622}, //    9 2019  1 10
  {58494.00, -0.0410966}, //   10 2019  1 11
  {58495.00, -0.0417283}, //   11 2019  1 12
  {58496.00, -0.0424595}, //   12 2019  1 13
  {58497.00, -0.0432639}, //   13 2019  1 14
  {58498.00, -0.0441369}, //   14 2019  1 15
  {58499.00, -0.0450440}, //   15 2019  1 16
  {58500.00, -0.0459161}, //   16 2019  1 17
  {58501.00, -0.0467200}, //   17 2019  1 18
  {58502.00, -0.0474747}, //   18 2019  1 19
  {58503.00, -0.0482159}, //   19 2019  1 20
  {58504.00, -0.0490102}, //   20 2019  1 21
  {58505.00, -0.0499233}, //   21 2019  1 22
  {58506.00, -0.0509608}, //   22 2019  1 23
  {58507.00, -0.0521836}, //   23 2019  1 24
  {58508.00, -0.0535609}, //   24 2019  1 25
  {58509.00, -0.0549981}, //   25 2019  1 26
  {58510.00, -0.0563823}, //   26 2019  1 27
  {58511.00, -0.0576299}, //   27 2019  1 28
  {58512.00, -0.0586946}, //   28 2019  1 29
  {58513.00, -0.0595854}, //   29 2019  1 30
  {58514.00, -0.0603394}, //   30 2019  1 31
  {58515.00, -0.0609865}, //   31 2019  2  1
  {58516.00, -0.0615645}, //   32 2019  2  2
  {58517.00, -0.0621140}, //   33 2019  2  3
  {58518.00, -0.0627000}, //   34 2019  2  4
  {58519.00, -0.0633608}, //   35 2019  2  5
  {58520.00, -0.0641101}, //   36 2019  2  6
  {58521.00, -0.0649794}, //   37 2019  2  7
  {58522.00, -0.0659661}, //   38 2019  2  8
  {58523.00, -0.0670379}, //   39 2019  2  9
  {58524.00, -0.0681559}, //   40 2019  2 10
  {58525.00, -0.0692735}, //   41 2019  2 11
  {58526.00, -0.0703556}, //   42 2019  2 12
  {58527.00, -0.0713787}, //   43 2019  2 13
  {58528.00, -0.0723210}, //   44 2019  2 14
  {58529.00, -0.0731898}, //   45 2019  2 15
  {58530.00, -0.0740111}, //   46 2019  2 16
  {58531.00, -0.0748515}, //   47 2019  2 17
  {58532.00, -0.0757862}, //   48 2019  2 18
  {58533.00, -0.0768900}, //   49 2019  2 19
  {58534.00, -0.0781927}, //   50 2019  2 20
  {58535.00, -0.0796674}, //   51 2019  2 21
  {58536.00, -0.0812286}, //   52 2019  2 22
  {58537.00, -0.0827714}, //   53 2019  2 23
  {58538.00, -0.0842068}, //   54 2019  2 24
  {58539.00, -0.0854693}, //   55 2019  2 25
  {58540.00, -0.0865362}, //   56 2019  2 26
  {58541.00, -0.0874277}, //   57 2019  2 27
  {58542.00, -0.0881818}, //   58 2019  2 28
  {58543.00, -0.0888453}, //   59 2019  3  1
  {58544.00, -0.0894670}, //   60 2019  3  2
  {58545.00, -0.0900991}, //   61 2019  3  3
  {58546.00, -0.0907970}, //   62 2019  3  4
  {58547.00, -0.0916016}, //   63 2019  3  5
  {58548.00, -0.0925374}, //   64 2019  3  6
  {58549.00, -0.0936081}, //   65 2019  3  7
  {58550.00, -0.0948003}, //   66 2019  3  8
  {58551.00, -0.0960791}, //   67 2019  3  9
  {58552.00, -0.0973977}, //   68 2019  3 10
  {58553.00, -0.0987052}, //   69 2019  3 11
  {58554.00, -0.0999588}, //   70 2019  3 12
  {58555.00, -0.1011315}, //   71 2019  3 13
  {58556.00, -0.1022138}, //   72 2019  3 14
  {58557.00, -0.1032225}, //   73 2019  3 15
  {58558.00, -0.1042053}, //   74 2019  3 16
  {58559.00, -0.1052390}, //   75 2019  3 17
  {58560.00, -0.1064108}, //   76 2019  3 18
  {58561.00, -0.1077893}, //   77 2019  3 19
  {58562.00, -0.1093868}, //   78 2019  3 20
  {58563.00, -0.1111512}, //   79 2019  3 21
  {58564.00, -0.1129795}, //   80 2019  3 22
  {58565.00, -0.1147512}, //   81 2019  3 23
  {58566.00, -0.1163666}, //   82 2019  3 24
  {58567.00, -0.1177707}, //   83 2019  3 25
  {58568.00, -0.1189584}, //   84 2019  3 26
  {58569.00, -0.1199627}, //   85 2019  3 27
  {58570.00, -0.1208381}, //   86 2019  3 28
  {58571.00, -0.1216464}, //   87 2019  3 29
  {58572.00, -0.1224458}, //   88 2019  3 30
  {58573.00, -0.1232880}, //   89 2019  3 31
  {58574.00, -0.1242133}, //   90 2019  4  1
  {58575.00, -0.1252478}, //   91 2019  4  2
  {58576.00, -0.1263987}, //   92 2019  4  3
  {58577.00, -0.1276577}, //   93 2019  4  4
  {58578.00, -0.1289972}, //   94 2019  4  5
  {58579.00, -0.1303726}, //   95 2019  4  6
  {58580.00, -0.1317309}, //   96 2019  4  7
  {58581.00, -0.1330196}, //   97 2019  4  8
  {58582.00, -0.1342005}, //   98 2019  4  9
  {58583.00, -0.1352587}, //   99 2019  4 10
  {58584.00, -0.1362110}, //  100 2019  4 11
  {58585.00, -0.1371036}, //  101 2019  4 12
  {58586.00, -0.1380086}, //  102 2019  4 13
  {58587.00, -0.1390094}, //  103 2019  4 14
  {58588.00, -0.1401761}, //  104 2019  4 15
  {58589.00, -0.1415433}, //  105 2019  4 16
  {58590.00, -0.1430925}, //  106 2019  4 17
  {58591.00, -0.1447503}, //  107 2019  4 18
  {58592.00, -0.1464085}, //  108 2019  4 19
  {58593.00, -0.1479566}, //  109 2019  4 20
  {58594.00, -0.1493155}, //  110 2019  4 21
  {58595.00, -0.1504528}, //  111 2019  4 22
  {58596.00, -0.1513812}, //  112 2019  4 23
  {58597.00, -0.1521470}, //  113 2019  4 24
  {58598.00, -0.1528094}, //  114 2019  4 25
  {58599.00, -0.1534293}, //  115 2019  4 26
  {58600.00, -0.1540632}, //  116 2019  4 27
  {58601.00, -0.1547572}, //  117 2019  4 28
  {58602.00, -0.1555450}, //  118 2019  4 29
  {58603.00, -0.1564409}, //  119 2019  4 30
  {58604.00, -0.1574408}, //  120 2019  5  1
  {58605.00, -0.1585250}, //  121 2019  5  2
  {58606.00, -0.1596578}, //  122 2019  5  3
  {58607.00, -0.1607911}, //  123 2019  5  4
  {58608.00, -0.1618692}, //  124 2019  5  5
  {58609.00, -0.1628441}, //  125 2019  5  6
  {58610.00, -0.1636905}, //  126 2019  5  7
  {58611.00, -0.1644132}, //  127 2019  5  8
  {58612.00, -0.1650509}, //  128 2019  5  9
  {58613.00, -0.1656698}, //  129 2019  5 10
  {58614.00, -0.1663483}, //  130 2019  5 11
  {58615.00, -0.1671549}, //  131 2019  5 12
  {58616.00, -0.1681103}, //  132 2019  5 13
  {58617.00, -0.1692313}, //  133 2019  5 14
  {58618.00, -0.1704719}, //  134 2019  5 15
  {58619.00, -0.1717479}, //  135 2019  5 16
  {58620.00, -0.1729587}, //  136 2019  5 17
  {58621.00, -0.1740150}, //  137 2019  5 18
  {58622.00, -0.1748623}, //  138 2019  5 19
  {58623.00, -0.1754906}, //  139 2019  5 20
  {58624.00, -0.1759282}, //  140 2019  5 21
  {58625.00, -0.1762273}, //  141 2019  5 22
  {58626.00, -0.1764490}, //  142 2019  5 23
  {58627.00, -0.1766521}, //  143 2019  5 24
  {58628.00, -0.1768876}, //  144 2019  5 25
  {58629.00, -0.1771933}, //  145 2019  5 26
  {58630.00, -0.1775901}, //  146 2019  5 27
  {58631.00, -0.1780811}, //  147 2019  5 28
  {58632.00, -0.1786526}, //  148 2019  5 29
  {58633.00, -0.1792767}, //  149 2019  5 30
  {58634.00, -0.1799128}, //  150 2019  5 31
  {58635.00, -0.1805109}, //  151 2019  6  1
  {58636.00, -0.1810193}, //  152 2019  6  2
  {58637.00, -0.1813977}, //  153 2019  6  3
  {58638.00, -0.1816350}, //  154 2019  6  4
  {58639.00, -0.1817584}, //  155 2019  6  5
  {58640.00, -0.1818320}, //  156 2019  6  6
  {58641.00, -0.1819376}, //  157 2019  6  7
  {58642.00, -0.1821495}, //  158 2019  6  8
  {58643.00, -0.1825122}, //  159 2019  6  9
  {58644.00, -0.1830280}, //  160 2019  6 10
  {58645.00, -0.1836579}, //  161 2019  6 11
  {58646.00, -0.1843300}, //  162 2019  6 12
  {58647.00, -0.1849579}, //  163 2019  6 13
  {58648.00, -0.1854609}, //  164 2019  6 14
  {58649.00, -0.1857825}, //  165 2019  6 15
  {58650.00, -0.1859017}, //  166 2019  6 16
  {58651.00, -0.1858335}, //  167 2019  6 17
  {58652.00, -0.1856211}, //  168 2019  6 18
  {58653.00, -0.1853209}, //  169 2019  6 19
  {58654.00, -0.1849924}, //  170 2019  6 20
  {58655.00, -0.1846891}, //  171 2019  6 21
  {58656.00, -0.1844533}, //  172 2019  6 22
  {58657.00, -0.1843120}, //  173 2019  6 23
  {58658.00, -0.1842732}, //  174 2019  6 24
  {58659.00, -0.1843271}, //  175 2019  6 25
  {58660.00, -0.1844492}, //  176 2019  6 26
  {58661.00, -0.1846049}, //  177 2019  6 27
  {58662.00, -0.1847526}, //  178 2019  6 28
  {58663.00, -0.1848465}, //  179 2019  6 29
  {58664.00, -0.1848447}, //  180 2019  6 30
  {58665.00, -0.1847224}, //  181 2019  7  1
  {58666.00, -0.1844879}, //  182 2019  7  2
  {58667.00, -0.1841911}, //  183 2019  7  3
  {58668.00, -0.1839139}, //  184 2019  7  4
  {58669.00, -0.1837437}, //  185 2019  7  5
  {58670.00, -0.1837407}, //  186 2019  7  6
  {58671.00, -0.1839158}, //  187 2019  7  7
  {58672.00, -0.1842273}, //  188 2019  7  8
  {58673.00, -0.1845973}, //  189 2019  7  9
  {58674.00, -0.1849342}, //  190 2019  7 10
  {58675.00, -0.1851569}, //  191 2019  7 11
  {58676.00, -0.1852099}, //  192 2019  7 12
  {58677.00, -0.1850719}, //  193 2019  7 13
  {58678.00, -0.1847543}, //  194 2019  7 14
  {58679.00, -0.1842943}, //  195 2019  7 15
  {58680.00, -0.1837441}, //  196 2019  7 16
  {58681.00, -0.1831606}, //  197 2019  7 17
  {58682.00, -0.1825970}, //  198 2019  7 18
  {58683.00, -0.1820968}, //  199 2019  7 19
  {58684.00, -0.1816897}, //  200 2019  7 20
  {58685.00, -0.1813879}, //  201 2019  7 21
  {58686.00, -0.1811849}, //  202 2019  7 22
  {58687.00, -0.1810581}, //  203 2019  7 23
  {58688.00, -0.1809732}, //  204 2019  7 24
  {58689.00, -0.1808908}, //  205 2019  7 25
  {58690.00, -0.1807701}, //  206 2019  7 26
  {58691.00, -0.1805728}, //  207 2019  7 27
  {58692.00, -0.1802703}, //  208 2019  7 28
  {58693.00, -0.1798579}, //  209 2019  7 29
  {58694.00, -0.1793684}, //  210 2019  7 30
  {58695.00, -0.1788763}, //  211 2019  7 31
  {58696.00, -0.1784821}, //  212 2019  8  1
  {58697.00, -0.1782795}, //  213 2019  8  2
  {58698.00, -0.1783191}, //  214 2019  8  3
  {58699.00, -0.1785879}, //  215 2019  8  4
  {58700.00, -0.1790144}, //  216 2019  8  5
  {58701.00, -0.1794955}, //  217 2019  8  6
  {58702.00, -0.1799310}, //  218 2019  8  7
  {58703.00, -0.1802497}, //  219 2019  8  8
  {58704.00, -0.1804202}, //  220 2019  8  9
  {58705.00, -0.1804477}, //  221 2019  8 10
  {58706.00, -0.1803641}, //  222 2019  8 11
  {58707.00, -0.1802169}, //  223 2019  8 12
  {58708.00, -0.1800592}, //  224 2019  8 13
  {58709.00, -0.1799423}, //  225 2019  8 14
  {58710.00, -0.1799088}, //  226 2019  8 15
  {58711.00, -0.1799877}, //  227 2019  8 16
  {58712.00, -0.1801899}, //  228 2019  8 17
  {58713.00, -0.1805066}, //  229 2019  8 18
  {58714.00, -0.1809106}, //  230 2019  8 19
  {58715.00, -0.1813635}, //  231 2019  8 20
  {58716.00, -0.1818224}, //  232 2019  8 21
  {58717.00, -0.1822462}, //  233 2019  8 22
  {58718.00, -0.1825997}, //  234 2019  8 23
  {58719.00, -0.1828606}, //  235 2019  8 24
  {58720.00, -0.1830265}, //  236 2019  8 25
  {58721.00, -0.1831164}, //  237 2019  8 26
  {58722.00, -0.1831690}, //  238 2019  8 27
  {58723.00, -0.1832451}, //  239 2019  8 28
  {58724.00, -0.1834238}, //  240 2019  8 29
  {58725.00, -0.1837750}, //  241 2019  8 30
  {58726.00, -0.1843141}, //  242 2019  8 31
  {58727.00, -0.1849832}, //  243 2019  9  1
  {58728.00, -0.1856947}, //  244 2019  9  2
  {58729.00, -0.1863463}, //  245 2019  9  3
  {58730.00, -0.1868555}, //  246 2019  9  4
  {58731.00, -0.1871843}, //  247 2019  9  5
  {58732.00, -0.1873367}, //  248 2019  9  6
  {58733.00, -0.1873501}, //  249 2019  9  7
  {58734.00, -0.1871656}, //  250 2019  9  8
  {58735.00, -0.1869129}, //  251 2019  9  9
  {58736.00, -0.1867204}, //  252 2019  9 10
  {58737.00, -0.1865911}, //  253 2019  9 11
  {58738.00, -0.1866391}, //  254 2019  9 12
  {58739.00, -0.1868626}, //  255 2019  9 13
  {58740.00, -0.1872800}, //  256 2019  9 14
  {58741.00, -0.1878722}, //  257 2019  9 15
  {58742.00, -0.1885847}, //  258 2019  9 16
  {58743.00, -0.1893267}, //  259 2019  9 17
  {58744.00, -0.1900513}, //  260 2019  9 18
  {58745.00, -0.1907722}, //  261 2019  9 19
  {58746.00, -0.1914382}, //  262 2019  9 20
  {58747.00, -0.1920677}, //  263 2019  9 21
  {58748.00, -0.1926597}, //  264 2019  9 22
  {58749.00, -0.1932348}, //  265 2019  9 23
  {58750.00, -0.1938773}, //  266 2019  9 24
  {58751.00, -0.1946966}, //  267 2019  9 25
  {58752.00, -0.1957595}, //  268 2019  9 26
  {58753.00, -0.1970383}, //  269 2019  9 27
  {58754.00, -0.1985664}, //  270 2019  9 28
  {58755.00, -0.2002747}, //  271 2019  9 29
  {58756.00, -0.2020003}, //  272 2019  9 30
  {58757.00, -0.2036559}, //  273 2019 10  1
  {58758.00, -0.2051065}, //  274 2019 10  2
  {58759.00, -0.2063310}, //  275 2019 10  3
  {58760.00, -0.2073173}, //  276 2019 10  4
  {58761.00, -0.2081454}, //  277 2019 10  5
  {58762.00, -0.2088326}, //  278 2019 10  6
  {58763.00, -0.2094856}, //  279 2019 10  7
  {58764.00, -0.2102273}, //  280 2019 10  8
  {58765.00, -0.2110946}, //  281 2019 10  9
  {58766.00, -0.2120226}, //  282 2019 10 10
  {58767.00, -0.2130144}, //  283 2019 10 11
  {58768.00, -0.2141134}, //  284 2019 10 12
  {58769.00, -0.2152960}, //  285 2019 10 13
  {58770.00, -0.2164713}, //  286 2019 10 14
  {58771.00, -0.2175630}, //  287 2019 10 15
  {58772.00, -0.2185636}, //  288 2019 10 16
  {58773.00, -0.2194791}, //  289 2019 10 17
  {58774.00, -0.2202325}, //  290 2019 10 18
  {58775.00, -0.2208963}, //  291 2019 10 19
  {58776.00, -0.2214415}, //  292 2019 10 20
  {58777.00, -0.2219983}, //  293 2019 10 21
  {58778.00, -0.2226602}, //  294 2019 10 22
  {58779.00, -0.2235048}, //  295 2019 10 23
  {58780.00, -0.2246190}, //  296 2019 10 24
  {58781.00, -0.2259492}, //  297 2019 10 25
  {58782.00, -0.2275062}, //  298 2019 10 26
  {58783.00, -0.2292217}, //  299 2019 10 27
  {58784.00, -0.2309042}, //  300 2019 10 28
  {58785.00, -0.2324249}, //  301 2019 10 29
  {58786.00, -0.2337534}, //  302 2019 10 30
  {58787.00, -0.2348841}, //  303 2019 10 31
  {58788.00, -0.2358067}, //  304 2019 11  1
  {58789.00, -0.2365908}, //  305 2019 11  2
  {58790.00, -0.2373340}, //  306 2019 11  3
  {58791.00, -0.2380879}, //  307 2019 11  4
  {58792.00, -0.2388739}, //  308 2019 11  5
  {58793.00, -0.2397562}, //  309 2019 11  6
  {58794.00, -0.2408110}, //  310 2019 11  7
  {58795.00, -0.2419493}, //  311 2019 11  8
  {58796.00, -0.2431590}, //  312 2019 11  9
  {58797.00, -0.2444017}, //  313 2019 11 10
  {58798.00, -0.2456428}, //  314 2019 11 11
  {58799.00, -0.2467618}, //  315 2019 11 12
  {58800.00, -0.2477403}, //  316 2019 11 13
  {58801.00, -0.2485698}, //  317 2019 11 14
  {58802.00, -0.2492314}, //  318 2019 11 15
  {58803.00, -0.2497674}, //  319 2019 11 16
  {58804.00, -0.2502091}, //  320 2019 11 17
  {58805.00, -0.2507040}, //  321 2019 11 18
  {58806.00, -0.2513338}, //  322 2019 11 19
  {58807.00, -0.2520896}, //  323 2019 11 20
  {58808.00, -0.2529502}, //  324 2019 11 21
  {58809.00, -0.2539650}, //  325 2019 11 22
  {58810.00, -0.2551195}, //  326 2019 11 23
  {58811.00, -0.2562953}, //  327 2019 11 24
  {58812.00, -0.2573934}, //  328 2019 11 25
  {58813.00, -0.2583719}, //  329 2019 11 26
  {58814.00, -0.2591494}, //  330 2019 11 27
  {58815.00, -0.2597076}, //  331 2019 11 28
  {58816.00, -0.2600871}, //  332 2019 11 29
  {58817.00, -0.2603506}, //  333 2019 11 30
  {58818.00, -0.2605961}, //  334 2019 12  1
  {58819.00, -0.2609119}, //  335 2019 12  2
  {58820.00, -0.2613181}, //  336 2019 12  3
  {58821.00, -0.2618088}, //  337 2019 12  4
  {58822.00, -0.2623809}, //  338 2019 12  5
  {58823.00, -0.2630280}, //  339 2019 12  6
  {58824.00, -0.2636826}, //  340 2019 12  7
  {58825.00, -0.2643763}, //  341 2019 12  8
  {58826.00, -0.2650727}, //  342 2019 12  9
  {58827.00, -0.2656933}, //  343 2019 12 10
  {58828.00, -0.2661553}, //  344 2019 12 11
  {58829.00, -0.2665041}, //  345 2019 12 12
  {58830.00, -0.2667725}, //  346 2019 12 13
  {58831.00, -0.2670052}, //  347 2019 12 14
  {58832.00, -0.2672403}, //  348 2019 12 15
  {58833.00, -0.2675352}, //  349 2019 12 16
  {58834.00, -0.2679590}, //  350 2019 12 17
  {58835.00, -0.2685566}, //  351 2019 12 18
  {58836.00, -0.2693574}, //  352 2019 12 19
  {58837.00, -0.2703452}, //  353 2019 12 20
  {58838.00, -0.2714234}, //  354 2019 12 21
  {58839.00, -0.2724743}, //  355 2019 12 22
  {58840.00, -0.2734090}, //  356 2019 12 23
  {58841.00, -0.2741483}, //  357 2019 12 24
  {58842.00, -0.2746538}, //  358 2019 12 25
  {58843.00, -0.2749645}, //  359 2019 12 26
  {58844.00, -0.2751411}, //  360 2019 12 27
  {58845.00, -0.2752837}, //  361 2019 12 28
  {58846.00, -0.2754225}, //  362 2019 12 29
  {58847.00, -0.2755313}, //  363 2019 12 30
  {58848.00, -0.2757593}, //  364 2019 12 31
  {58849.00, -0.2761056}, //  365 2020  1  1
  {58850.00, -0.2765286}, //  366 2020  1  2
  {58851.00, -0.2770727}, //  367 2020  1  3
  {58852.00, -0.2776348}, //  368 2020  1  4
  {58853.00, -0.2781239}, //  369 2020  1  5
  {58854.00, -0.2785952}, //  370 2020  1  6
  {58855.00, -0.2790371}, //  371 2020  1  7
  {58856.00, -0.2793816}, //  372 2020  1  8
  {58857.00, -0.2796173}, //  373 2020  1  9
  {58858.00, -0.2798010}, //  374 2020  1 10
  {58859.00, -0.2800178}, //  375 2020  1 11
  {58860.00, -0.2803257}, //  376 2020  1 12
  {58861.00, -0.2808219}, //  377 2020  1 13
  {58862.00, -0.2815751}, //  378 2020  1 14
  {58863.00, -0.2825685}, //  379 2020  1 15
  {58864.00, -0.2837300}, //  380 2020  1 16
  {58865.00, -0.2849780}, //  381 2020  1 17
  {58866.00, -0.2862136}, //  382 2020  1 18
  {58867.00, -0.2873665}, //  383 2020  1 19
  {58868.00, -0.2883673}, //  384 2020  1 20
  {58869.00, -0.2892033}, //  385 2020  1 21
  {58870.00, -0.2898792}, //  386 2020  1 22
  {58871.00, -0.2904304}, //  387 2020  1 23
  {58872.00, -0.2909817}, //  388 2020  1 24
  {58873.00, -0.2916560}, //  389 2020  1 25
  {58874.00, -0.2924478}, //  390 2020  1 26
  {58875.00, -0.2933674}, //  391 2020  1 27
  {58876.00, -0.2944534}, //  392 2020  1 28
  {58877.00, -0.2956681}, //  393 2020  1 29
  {58878.00, -0.2969503}, //  394 2020  1 30
  {58879.00, -0.2982275}, //  395 2020  1 31
  {58880.00, -0.2994122}  //  396 2020  2  1
};
#endif /* eopdut_h */