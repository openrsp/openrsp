/* OpenRSP: open-ended library for response theory
   Copyright 2014

   OpenRSP is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   OpenRSP is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with OpenRSP. If not, see <http://www.gnu.org/licenses/>.

   This file initializes unperturbed atomic orbital Fock matrix.

   2015-02-10, Bin Gao:
   * first version
*/

#include "tests/ao_dens/openrsp_c_ao_fock.h"

const QReal values_fock[] = {
      -20.632274495919347, -5.141049902525445, -3.649348247710690,
       -0.007051562069756,  0.014889563995715, -0.018454432509332,
       -0.001487288582570,  0.002659657423952, -0.003899729961522,
        0.000394876038588,  0.000812344307960, -0.001295838511869,
        0.001327337799432, -0.004967167815125, -0.000078231317327,
       -0.010789705296783, -0.466673088551037,  0.000041514574558,
       -0.031579702367395, -0.000031186469981, -0.000249017518649,
       -1.311114096489859, -0.000605504077390, -0.000108369486437,
       -0.000439474777299,  0.148244855027145,  0.000010664620884,
        0.256643897818593, -0.742781038795938, -1.496861365782057,
       -0.002735187363306, -0.324204582428558, -5.141049902525445,
       -2.167965660236860, -2.129010392187726, -0.015609587406240,
       -0.012298049230703, -0.043437569503949, -0.018470265771879,
        0.000805326126169, -0.050877522104937, -0.010358886170551,
       -0.019683262947666, -0.042741042573198, -0.045067396438292,
        0.021267402105457, -0.010789705296783, -0.131713636219224,
       -0.458761574680208, -0.001088063595550, -0.227055215610765,
       -0.001260565776554, -0.002224484000979, -0.961489322358230,
       -0.005264873760155,  0.001461435645005, -0.002054550987210,
        0.183071660341150, -0.000871448527813,  0.316633845043071,
       -0.795053128631225, -1.146364048792113, -0.028104370519575,
       -0.317262826959862, -3.649348247710690, -2.129010392187726,
       -1.706347875174290, -0.028484057909972,  0.021012904565641,
       -0.077130737418269, -0.050927899845191,  0.162118970323373,
       -0.148200335053003, -0.017945640405865, -0.031120349268672,
       -0.046748880095479, -0.077041246604314,  0.083045495601476,
       -0.466673088551036, -0.458761574680205, -0.632800990619478,
        0.003792772371608, -0.278989596073572, -0.012139410181443,
        0.000505082898453, -0.754567806756491, -0.038836261148490,
        0.000823762456693, -0.013962600309319,  0.098108152716691,
        0.007169374298154,  0.189695730308029, -0.709286545234027,
       -0.977341379441460, -0.146895315564143, -0.414634223828795,
       -0.007051562069756, -0.015609587406240, -0.028484057909972,
       -0.241869441973902,  0.007089388012289,  0.021560715739354,
       -0.680124398116337, -0.000126067409724, -0.010500536157708,
        0.013935204103084, -0.015260569859204, -0.025569875758014,
       -0.087831988899201, -0.025877109558907, -0.000041514574558,
        0.001088063595549, -0.003792772371608, -0.034991133076481,
       -0.000753429375744,  0.000512623793870, -0.175151248269819,
       -0.010914132693093,  0.000783557108506, -0.148565346712821,
        0.005384886135021, -0.001361773627893, -0.000739612026537,
       -0.000102456067582, -0.237364783628079, -0.137431881689734,
        0.005777294935823,  0.025768972174841,  0.014889563995715,
       -0.012298049230703,  0.021012904565641,  0.007089388012289,
       -0.005189999477224,  0.015156872831067, -0.000231667854541,
       -0.598492635548030,  0.001555594441078, -0.018830932196900,
       -0.059365797610609, -0.081074871317186, -0.013490172993269,
       -0.115164583304221,  0.031579702367395,  0.227055215610768,
        0.278989596073573, -0.000753429375743,  0.357737841551756,
        0.001214448755987,  0.002307518132586,  0.349788405546355,
        0.004867550650933, -0.000132181557939,  0.014810420506377,
       -0.104919405793600,  0.002341671257552, -0.177780274814989,
       -0.098477074749993, -0.035213719777713,  0.037303797059004,
        0.135018163484147, -0.018454432509332, -0.043437569503949,
       -0.077130737418269,  0.021560715739354,  0.015156872831067,
       -0.199720413951551, -0.010622792727178,  0.001375798876128,
       -0.704406115882429, -0.014604965083661, -0.010680648240933,
       -0.154645686726750, -0.091533407711460, -0.010825180236378,
       -0.000031186469981, -0.001260565776554, -0.012139410181443,
       -0.000512623793870, -0.001214448755988, -0.037754894990004,
        0.000139739603144, -0.031109789884863, -0.174971978432339,
        0.004912885572421, -0.142162057075709,  0.000070952545540,
       -0.003112358516166,  0.002603202794562, -0.573217834870032,
       -0.336633710144226, -0.017149185135695, -0.088515593771827,
       -0.001487288582570, -0.018470265771879, -0.050927899845191,
       -0.680124398116337, -0.000231667854541, -0.010622792727178,
       -0.037258194907669, -0.015017299904310, -0.057135441655455,
        0.032184234705792, -0.020244225254748, -0.039258891526505,
       -0.106386091252830, -0.029903618075322,  0.000249017518649,
        0.002224484000979, -0.000505082898453, -0.175151248269817,
        0.002307518132586, -0.000139739603144, -0.117360734477403,
       -0.015043483288293,  0.000249725589133, -0.014354515964218,
        0.001362959663099,  0.004023549318176, -0.016629595038010,
        0.004130142200289, -0.162238910357895, -0.090683929477314,
        0.044417127363252,  0.019216280917415,  0.002659657423952,
        0.000805326126169,  0.162118970323373, -0.000126067409724,
       -0.598492635548030,  0.001375798876128, -0.015017299904310,
       -0.333014834611723, -0.012185202279149, -0.021722268691638,
       -0.074824524080580, -0.147925103682673, -0.020685868772494,
       -0.220243710186532,  1.311114096489854,  0.961489322358220,
        0.754567806756486, -0.010914132693093,  0.349788405546359,
        0.031109789884864, -0.015043483288292,  0.113302575229062,
        0.058059178270948, -0.000358538725781,  0.022776912693907,
        0.020881524677648, -0.021482136301287, -0.013796940078492,
       -0.064791173891956,  0.091140202485237,  0.297628309631096,
        0.406729527296176, -0.003899729961522, -0.050877522104937,
       -0.148200335053003, -0.010500536157708,  0.001555594441078,
       -0.704406115882429, -0.057135441655455, -0.012185202279149,
       -0.166163383238192, -0.019932197579250, -0.002920862124709,
       -0.198455542065005, -0.120264911811128, -0.015490248106323,
       -0.000605504077390, -0.005264873760155, -0.038836261148490,
       -0.000783557108505, -0.004867550650933, -0.174971978432337,
       -0.000249725589133, -0.058059178270948, -0.145544676949275,
        0.001144385747621, -0.015029657631162, -0.026466609133454,
        0.012411439921194,  0.002833644530990, -0.397512788571889,
       -0.246308898397746, -0.110930883255261, -0.111773981148103,
        0.000394876038588, -0.010358886170551, -0.017945640405865,
        0.013935204103084, -0.018830932196900, -0.014604965083661,
        0.032184234705792, -0.021722268691638, -0.019932197579250,
        2.009583138274905, -0.024503807372761, -0.006659922473222,
       -0.017286991675712, -0.001846976952510, -0.000108369486437,
        0.001461435645005,  0.000823762456693,  0.148565346712821,
        0.000132181557940,  0.004912885572420,  0.014354515964219,
        0.000358538725782,  0.001144385747622, -0.108481500023493,
        0.000042180354717, -0.002146635635280,  0.007656397119497,
       -0.002307176838065,  0.002540964315591, -0.014822004710886,
       -0.029488700950358, -0.009830449683270,  0.000812344307960,
       -0.019683262947666, -0.031120349268672, -0.015260569859204,
       -0.059365797610609, -0.010680648240933, -0.020244225254748,
       -0.074824524080580, -0.002920862124709, -0.024503807372761,
        1.956468166118409, -0.034828951110692, -0.029003985388192,
        0.000820456906258,  0.000439474777299,  0.002054550987211,
        0.013962600309320,  0.005384886135021,  0.014810420506378,
        0.142162057075710,  0.001362959663099,  0.022776912693907,
        0.015029657631163, -0.000042180354717, -0.106011526341941,
        0.005411120187539, -0.004863647378398, -0.011510713839348,
        0.009040559553115, -0.027244823845372,  0.073578854480256,
        0.028062514420651, -0.001295838511869, -0.042741042573198,
       -0.046748880095479, -0.025569875758014, -0.081074871317186,
       -0.154645686726750, -0.039258891526505, -0.147925103682673,
       -0.198455542065005, -0.006659922473222, -0.034828951110692,
        1.881713044271649, -0.075483879749266, -0.082008409369252,
        0.148244855027145,  0.183071660341150,  0.098108152716691,
        0.001361773627893,  0.104919405793601,  0.000070952545540,
       -0.004023549318176, -0.020881524677646, -0.026466609133454,
       -0.002146635635280, -0.005411120187539,  0.090166018360330,
       -0.001002459533756,  0.177339331130028, -0.011488478691420,
       -0.105835638851651,  0.034355519785872,  0.035248189544819,
        0.001327337799432, -0.045067396438292, -0.077041246604314,
       -0.087831988899201, -0.013490172993269, -0.091533407711460,
       -0.106386091252830, -0.020685868772494, -0.120264911811128,
       -0.017286991675712, -0.029003985388192, -0.075483879749266,
        1.966474207685003, -0.029489149756268, -0.000010664620884,
        0.000871448527813, -0.007169374298154, -0.000739612026537,
        0.002341671257552,  0.003112358516166, -0.016629595038010,
       -0.021482136301287, -0.012411439921194, -0.007656397119497,
       -0.004863647378398,  0.001002459533756, -0.013754999455413,
        0.001261249932155, -0.004450381386899, -0.102040867867434,
        0.010674020208523, -0.013106057332217, -0.004967167815125,
        0.021267402105457,  0.083045495601476, -0.025877109558907,
       -0.115164583304221, -0.010825180236378, -0.029903618075322,
       -0.220243710186532, -0.015490248106323, -0.001846976952510,
        0.000820456906258, -0.082008409369252, -0.029489149756268,
        1.936996498941165,  0.256643897818593,  0.316633845043074,
        0.189695730308025,  0.000102456067582,  0.177780274814998,
        0.002603202794562, -0.004130142200289,  0.013796940078491,
        0.002833644530990, -0.002307176838065,  0.011510713839348,
        0.177339331130029, -0.001261249932154,  0.290335647894065,
       -0.002979670330502,  0.032493380318085,  0.080375332494570,
        0.106249857795518, -0.000078231317327, -0.010789705296783,
       -0.466673088551036, -0.000041514574558,  0.031579702367395,
       -0.000031186469981,  0.000249017518649,  1.311114096489854,
       -0.000605504077390, -0.000108369486437,  0.000439474777299,
        0.148244855027145, -0.000010664620884,  0.256643897818593,
      -20.632274495919241, -5.141049902525430, -3.649348247710672,
        0.007051562069755, -0.014889563995715, -0.018454432509332,
        0.001487288582570, -0.002659657423952, -0.003899729961522,
        0.000394876038588, -0.000812344307960, -0.001295838511868,
       -0.001327337799432, -0.004967167815124, -0.002735187363306,
       -0.324204582428557, -0.742781038795936, -1.496861365782052,
       -0.010789705296783, -0.131713636219224, -0.458761574680205,
        0.001088063595549,  0.227055215610768, -0.001260565776554,
        0.002224484000979,  0.961489322358220, -0.005264873760155,
        0.001461435645005,  0.002054550987211,  0.183071660341150,
        0.000871448527813,  0.316633845043074, -5.141049902525430,
       -2.167965660236822, -2.129010392187689,  0.015609587406241,
        0.012298049230695, -0.043437569503950,  0.018470265771880,
       -0.000805326126176, -0.050877522104937, -0.010358886170551,
        0.019683262947666, -0.042741042573194,  0.045067396438292,
        0.021267402105461, -0.028104370519575, -0.317262826959857,
       -0.795053128631216, -1.146364048792092, -0.466673088551037,
       -0.458761574680208, -0.632800990619478, -0.003792772371608,
        0.278989596073573, -0.012139410181443, -0.000505082898453,
        0.754567806756486, -0.038836261148490,  0.000823762456693,
        0.013962600309320,  0.098108152716691, -0.007169374298154,
        0.189695730308025, -3.649348247710672, -2.129010392187689,
       -1.706347875174240,  0.028484057909973, -0.021012904565647,
       -0.077130737418270,  0.050927899845191, -0.162118970323383,
       -0.148200335053003, -0.017945640405866,  0.031120349268672,
       -0.046748880095477,  0.077041246604314,  0.083045495601480,
       -0.146895315564143, -0.414634223828791, -0.709286545234017,
       -0.977341379441445,  0.000041514574558, -0.001088063595550,
        0.003792772371608, -0.034991133076481, -0.000753429375743,
       -0.000512623793870, -0.175151248269817, -0.010914132693093,
       -0.000783557108505,  0.148565346712821,  0.005384886135021,
        0.001361773627893, -0.000739612026537,  0.000102456067582,
        0.007051562069755,  0.015609587406241,  0.028484057909973,
       -0.241869441973858,  0.007089388012290, -0.021560715739350,
       -0.680124398116335, -0.000126067409723,  0.010500536157709,
       -0.013935204103088, -0.015260569859203,  0.025569875758014,
       -0.087831988899200,  0.025877109558907, -0.005777294935823,
       -0.025768972174840,  0.237364783628074,  0.137431881689731,
       -0.031579702367395, -0.227055215610765, -0.278989596073572,
       -0.000753429375744,  0.357737841551756, -0.001214448755988,
        0.002307518132586,  0.349788405546359, -0.004867550650933,
        0.000132181557940,  0.014810420506378,  0.104919405793601,
        0.002341671257552,  0.177780274814998, -0.014889563995715,
        0.012298049230695, -0.021012904565647,  0.007089388012290,
       -0.005189999477242, -0.015156872831067, -0.000231667854541,
       -0.598492635548026, -0.001555594441078,  0.018830932196899,
       -0.059365797610610,  0.081074871317189, -0.013490172993270,
        0.115164583304226, -0.037303797059004, -0.135018163484148,
        0.098477074749992,  0.035213719777710, -0.000031186469981,
       -0.001260565776554, -0.012139410181443,  0.000512623793870,
        0.001214448755987, -0.037754894990004, -0.000139739603144,
        0.031109789884864, -0.174971978432337,  0.004912885572420,
        0.142162057075710,  0.000070952545540,  0.003112358516166,
        0.002603202794562, -0.018454432509332, -0.043437569503950,
       -0.077130737418270, -0.021560715739350, -0.015156872831067,
       -0.199720413951463,  0.010622792727179, -0.001375798876128,
       -0.704406115882415, -0.014604965083661,  0.010680648240928,
       -0.154645686726750,  0.091533407711459, -0.010825180236378,
       -0.017149185135695, -0.088515593771826, -0.573217834870019,
       -0.336633710144224, -0.000249017518649, -0.002224484000979,
        0.000505082898453, -0.175151248269819,  0.002307518132586,
        0.000139739603144, -0.117360734477403, -0.015043483288292,
       -0.000249725589133,  0.014354515964219,  0.001362959663099,
       -0.004023549318176, -0.016629595038010, -0.004130142200289,
        0.001487288582570,  0.018470265771880,  0.050927899845191,
       -0.680124398116335, -0.000231667854541,  0.010622792727179,
       -0.037258194907682, -0.015017299904310,  0.057135441655455,
       -0.032184234705799, -0.020244225254748,  0.039258891526505,
       -0.106386091252829,  0.029903618075322, -0.044417127363253,
       -0.019216280917414,  0.162238910357886,  0.090683929477309,
       -1.311114096489859, -0.961489322358230, -0.754567806756491,
       -0.010914132693093,  0.349788405546355, -0.031109789884863,
       -0.015043483288293,  0.113302575229062, -0.058059178270948,
        0.000358538725782,  0.022776912693907, -0.020881524677646,
       -0.021482136301287,  0.013796940078491, -0.002659657423952,
       -0.000805326126176, -0.162118970323383, -0.000126067409723,
       -0.598492635548026, -0.001375798876128, -0.015017299904310,
       -0.333014834611713,  0.012185202279149,  0.021722268691638,
       -0.074824524080580,  0.147925103682676, -0.020685868772494,
        0.220243710186537, -0.297628309631099, -0.406729527296191,
        0.064791173891951, -0.091140202485245, -0.000605504077390,
       -0.005264873760155, -0.038836261148490,  0.000783557108506,
        0.004867550650933, -0.174971978432339,  0.000249725589133,
        0.058059178270948, -0.145544676949275,  0.001144385747622,
        0.015029657631163, -0.026466609133454, -0.012411439921194,
        0.002833644530990, -0.003899729961522, -0.050877522104937,
       -0.148200335053003,  0.010500536157709, -0.001555594441078,
       -0.704406115882415,  0.057135441655455,  0.012185202279149,
       -0.166163383238162, -0.019932197579250,  0.002920862124702,
       -0.198455542065005,  0.120264911811128, -0.015490248106323,
       -0.110930883255262, -0.111773981148100, -0.397512788571874,
       -0.246308898397741, -0.000108369486437,  0.001461435645005,
        0.000823762456693, -0.148565346712821, -0.000132181557939,
        0.004912885572421, -0.014354515964218, -0.000358538725781,
        0.001144385747621, -0.108481500023493, -0.000042180354717,
       -0.002146635635280, -0.007656397119497, -0.002307176838065,
        0.000394876038588, -0.010358886170551, -0.017945640405866,
       -0.013935204103088,  0.018830932196899, -0.014604965083661,
       -0.032184234705799,  0.021722268691638, -0.019932197579250,
        2.009583138274929,  0.024503807372761, -0.006659922473222,
        0.017286991675712, -0.001846976952510, -0.029488700950358,
       -0.009830449683270,  0.002540964315593, -0.014822004710885,
       -0.000439474777299, -0.002054550987210, -0.013962600309319,
        0.005384886135021,  0.014810420506377, -0.142162057075709,
        0.001362959663099,  0.022776912693907, -0.015029657631162,
        0.000042180354717, -0.106011526341941, -0.005411120187539,
       -0.004863647378398,  0.011510713839348, -0.000812344307960,
        0.019683262947666,  0.031120349268672, -0.015260569859203,
       -0.059365797610610,  0.010680648240928, -0.020244225254748,
       -0.074824524080580,  0.002920862124702,  0.024503807372761,
        1.956468166118430,  0.034828951110692, -0.029003985388192,
       -0.000820456906258, -0.073578854480256, -0.028062514420650,
       -0.009040559553123,  0.027244823845369,  0.148244855027145,
        0.183071660341150,  0.098108152716691, -0.001361773627893,
       -0.104919405793600,  0.000070952545540,  0.004023549318176,
        0.020881524677648, -0.026466609133454, -0.002146635635280,
        0.005411120187539,  0.090166018360330,  0.001002459533756,
        0.177339331130029, -0.001295838511868, -0.042741042573194,
       -0.046748880095477,  0.025569875758014,  0.081074871317189,
       -0.154645686726750,  0.039258891526505,  0.147925103682676,
       -0.198455542065005, -0.006659922473222,  0.034828951110692,
        1.881713044271702,  0.075483879749266, -0.082008409369252,
        0.034355519785872,  0.035248189544820, -0.011488478691422,
       -0.105835638851648,  0.000010664620884, -0.000871448527813,
        0.007169374298154, -0.000739612026537,  0.002341671257552,
       -0.003112358516166, -0.016629595038010, -0.021482136301287,
        0.012411439921194,  0.007656397119497, -0.004863647378398,
       -0.001002459533756, -0.013754999455413, -0.001261249932154,
       -0.001327337799432,  0.045067396438292,  0.077041246604314,
       -0.087831988899200, -0.013490172993270,  0.091533407711459,
       -0.106386091252829, -0.020685868772494,  0.120264911811128,
        0.017286991675712, -0.029003985388192,  0.075483879749266,
        1.966474207685044,  0.029489149756268, -0.010674020208523,
        0.013106057332217,  0.004450381386886,  0.102040867867432,
        0.256643897818593,  0.316633845043071,  0.189695730308029,
       -0.000102456067582, -0.177780274814989,  0.002603202794562,
        0.004130142200289, -0.013796940078492,  0.002833644530990,
       -0.002307176838065, -0.011510713839348,  0.177339331130028,
        0.001261249932155,  0.290335647894065, -0.004967167815124,
        0.021267402105461,  0.083045495601480,  0.025877109558907,
        0.115164583304226, -0.010825180236378,  0.029903618075322,
        0.220243710186537, -0.015490248106323, -0.001846976952510,
       -0.000820456906258, -0.082008409369252,  0.029489149756268,
        1.936996498941182,  0.080375332494571,  0.106249857795519,
       -0.002979670330500,  0.032493380318088, -0.742781038795938,
       -0.795053128631225, -0.709286545234027, -0.237364783628079,
       -0.098477074749993, -0.573217834870032, -0.162238910357895,
       -0.064791173891956, -0.397512788571889,  0.002540964315591,
        0.009040559553115, -0.011488478691420, -0.004450381386899,
       -0.002979670330502, -0.002735187363306, -0.028104370519575,
       -0.146895315564143, -0.005777294935823, -0.037303797059004,
       -0.017149185135695, -0.044417127363253, -0.297628309631099,
       -0.110930883255262, -0.029488700950358, -0.073578854480256,
        0.034355519785872, -0.010674020208523,  0.080375332494571,
       -0.254678826425118, -0.648515327506478, -0.025937685980667,
       -0.165330764433261, -1.496861365782057, -1.146364048792113,
       -0.977341379441460, -0.137431881689734, -0.035213719777713,
       -0.336633710144226, -0.090683929477314,  0.091140202485237,
       -0.246308898397746, -0.014822004710886, -0.027244823845372,
       -0.105835638851651, -0.102040867867434,  0.032493380318085,
       -0.324204582428557, -0.317262826959857, -0.414634223828791,
       -0.025768972174840, -0.135018163484148, -0.088515593771826,
       -0.019216280917414, -0.406729527296191, -0.111773981148100,
       -0.009830449683270, -0.028062514420650,  0.035248189544820,
        0.013106057332217,  0.106249857795519, -0.648515327506478,
       -0.618798048851664, -0.165330764433259, -0.310034259569418,
       -0.002735187363306, -0.028104370519575, -0.146895315564143,
        0.005777294935823,  0.037303797059004, -0.017149185135695,
        0.044417127363252,  0.297628309631096, -0.110930883255261,
       -0.029488700950358,  0.073578854480256,  0.034355519785872,
        0.010674020208523,  0.080375332494570, -0.742781038795936,
       -0.795053128631216, -0.709286545234017,  0.237364783628074,
        0.098477074749992, -0.573217834870019,  0.162238910357886,
        0.064791173891951, -0.397512788571874,  0.002540964315593,
       -0.009040559553123, -0.011488478691422,  0.004450381386886,
       -0.002979670330500, -0.025937685980667, -0.165330764433259,
       -0.254678826425089, -0.648515327506461, -0.324204582428558,
       -0.317262826959862, -0.414634223828795,  0.025768972174841,
        0.135018163484147, -0.088515593771827,  0.019216280917415,
        0.406729527296176, -0.111773981148103, -0.009830449683270,
        0.028062514420651,  0.035248189544819, -0.013106057332217,
        0.106249857795518, -1.496861365782052, -1.146364048792092,
       -0.977341379441445,  0.137431881689731,  0.035213719777710,
       -0.336633710144224,  0.090683929477309, -0.091140202485245,
       -0.246308898397741, -0.014822004710885,  0.027244823845369,
       -0.105835638851648,  0.102040867867432,  0.032493380318088,
       -0.165330764433261, -0.310034259569418, -0.648515327506461,
       -0.618798048851638
};
