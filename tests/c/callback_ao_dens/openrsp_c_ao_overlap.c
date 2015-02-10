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

   This file initializes overlap integrals.

   2015-02-10, Bin Gao:
   * first version
*/

#include "tests/ao_dens/openrsp_c_ao_overlap.h"

const QReal values_overlap[] = {
        1.000000000000000,  0.233689906059163,  0.167279759595895,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000013322,
        0.000360781306443,  0.021035335413722,  0.000000000000000,
        0.001155817773709,  0.000000000000000,  0.000000000000000,
        0.059605849715634,  0.000000000000000,  0.000000000000000,
        0.000000000000000, -0.006457144066568,  0.000000000000000,
       -0.011184101595088,  0.033253425203908,  0.068076201809780,
        0.000104734953463,  0.014626666566694,  0.233689906059163,
        1.000000000000000,  0.763640773338444,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000360781306443,  0.022466373968732,
        0.144341571008982,  0.000000000000000,  0.042457779708000,
        0.000000000000000,  0.000000000000000,  0.329505331313828,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
       -0.070665963942989,  0.000000000000000, -0.122397039915088,
        0.240404328848086,  0.375934821334508,  0.005168519068412,
        0.098291977680733,  0.167279759595895,  0.763640773338444,
        1.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.021035335413722,  0.144341571008982,  0.342899055855409,
        0.000000000000000,  0.130703157963586,  0.000000000000000,
        0.000000000000000,  0.501692882467059,  0.000000000000000,
        0.000000000000000,  0.000000000000000, -0.076197708983921,
        0.000000000000000, -0.131978303380498,  0.420483682453648,
        0.678093402858836,  0.063349060947094,  0.254121532334353,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        1.000000000000000,  0.000000000000000,  0.000000000000000,
        0.501520641074811,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.010097520861493,
        0.000000000000000,  0.000000000000000,  0.089333516637993,
        0.000000000000000,  0.000000000000000,  0.062253543338693,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.105654271868166,  0.055743215005449,
       -0.001871814694006, -0.014195841148136,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        1.000000000000000,  0.000000000000000,  0.000000000000000,
        0.501520641074811,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000, -0.001155817773709, -0.042457779708000,
       -0.130703157963586,  0.000000000000000, -0.073746808472246,
        0.000000000000000,  0.000000000000000, -0.218104134283380,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.095489194679018,  0.000000000000000,  0.165392136757895,
        0.048681748856291,  0.025684500449954, -0.008491793599620,
       -0.064401755894416,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        1.000000000000000,  0.000000000000000,  0.000000000000000,
        0.501520641074811,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.010097520861493,
        0.000000000000000,  0.000000000000000,  0.089333516637993,
        0.000000000000000,  0.062253543338693,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.255056682446983,
        0.134567956759564,  0.004518689472437,  0.034269737359089,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.501520641074811,  0.000000000000000,  0.000000000000000,
        1.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.089333516637993,
        0.000000000000000,  0.000000000000000,  0.342899055855409,
        0.000000000000000,  0.000000000000000,  0.180410116244944,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.219175793213426,  0.182045269436310,
       -0.032607081177624, -0.068223083469529,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.501520641074811,  0.000000000000000,  0.000000000000000,
        1.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000, -0.059605849715634, -0.329505331313828,
       -0.501692882467059,  0.000000000000000, -0.218104134283380,
        0.000000000000000,  0.000000000000000, -0.391123811866453,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.062544706021092,  0.000000000000000,  0.108330608572991,
        0.100988447811255,  0.083880016685293, -0.147927358478970,
       -0.309505179870655,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.501520641074811,  0.000000000000000,  0.000000000000000,
        1.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.089333516637993,
        0.000000000000000,  0.000000000000000,  0.342899055855409,
        0.000000000000000,  0.180410116244944,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.529105446483568,
        0.439469807103769,  0.078715737682824,  0.164695217981853,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        1.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000, -0.062253543338693,
        0.000000000000000,  0.000000000000000, -0.180410116244944,
        0.000000000000000,  0.000000000000000, -0.224124179944020,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.043677731061311,  0.008151460212183,
        0.010755059186520,  0.021770630621996,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        1.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
       -0.062253543338693,  0.000000000000000,  0.000000000000000,
       -0.180410116244944,  0.000000000000000, -0.224124179944020,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.105441048282554,  0.019678185860883, -0.025963452940715,
       -0.052555800376472,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        1.000000000000001,  0.000000000000000,  0.000000000000000,
       -0.006457144066568, -0.070665963942989, -0.076197708983921,
        0.000000000000000, -0.095489194679018,  0.000000000000000,
        0.000000000000000, -0.062544706021092,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.093891255171399,
        0.000000000000000,  0.089962792984985,  0.285773499802488,
        0.053333157577861, -0.006792877136544, -0.013750293367576,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        1.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.041951212424467,
        0.000000000000000,  0.228839296924101,  0.042707676853505,
       -0.005723028021282, -0.011584680932934,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        1.000000000000000, -0.011184101595088, -0.122397039915088,
       -0.131978303380498,  0.000000000000000, -0.165392136757896,
        0.000000000000000,  0.000000000000000, -0.108330608572991,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.089962792984985,  0.000000000000000,  0.197771340665263,
        0.037334435775198,  0.006967627675958, -0.023210693813002,
       -0.046983603969048,  0.000000000013322,  0.000360781306443,
        0.021035335413722,  0.000000000000000, -0.001155817773709,
        0.000000000000000,  0.000000000000000, -0.059605849715634,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
       -0.006457144066568,  0.000000000000000, -0.011184101595088,
        1.000000000000000,  0.233689906059163,  0.167279759595895,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000104734953463,
        0.014626666566694,  0.033253425203908,  0.068076201809780,
        0.000360781306443,  0.022466373968732,  0.144341571008982,
        0.000000000000000, -0.042457779708000,  0.000000000000000,
        0.000000000000000, -0.329505331313828,  0.000000000000000,
        0.000000000000000,  0.000000000000000, -0.070665963942989,
        0.000000000000000, -0.122397039915088,  0.233689906059163,
        1.000000000000000,  0.763640773338444,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.005168519068412,  0.098291977680733,
        0.240404328848086,  0.375934821334508,  0.021035335413722,
        0.144341571008982,  0.342899055855409,  0.000000000000000,
       -0.130703157963586,  0.000000000000000,  0.000000000000000,
       -0.501692882467059,  0.000000000000000,  0.000000000000000,
        0.000000000000000, -0.076197708983921,  0.000000000000000,
       -0.131978303380498,  0.167279759595895,  0.763640773338444,
        1.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.063349060947094,  0.254121532334353,  0.420483682453648,
        0.678093402858836,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.010097520861493,  0.000000000000000,
        0.000000000000000,  0.089333516637993,  0.000000000000000,
        0.000000000000000, -0.062253543338693,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        1.000000000000000,  0.000000000000000,  0.000000000000000,
        0.501520641074811,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.001871814694006,
        0.014195841148136, -0.105654271868166, -0.055743215005449,
        0.001155817773709,  0.042457779708000,  0.130703157963586,
        0.000000000000000, -0.073746808472246,  0.000000000000000,
        0.000000000000000, -0.218104134283380,  0.000000000000000,
        0.000000000000000,  0.000000000000000, -0.095489194679018,
        0.000000000000000, -0.165392136757896,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        1.000000000000000,  0.000000000000000,  0.000000000000000,
        0.501520641074811,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.008491793599620,  0.064401755894416,
       -0.048681748856291, -0.025684500449954,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.010097520861493,  0.000000000000000,
        0.000000000000000,  0.089333516637993,  0.000000000000000,
       -0.062253543338693,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        1.000000000000000,  0.000000000000000,  0.000000000000000,
        0.501520641074811,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.004518689472437,  0.034269737359089,  0.255056682446983,
        0.134567956759564,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.089333516637993,  0.000000000000000,
        0.000000000000000,  0.342899055855409,  0.000000000000000,
        0.000000000000000, -0.180410116244944,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.501520641074811,  0.000000000000000,  0.000000000000000,
        1.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.032607081177624,
        0.068223083469529, -0.219175793213426, -0.182045269436310,
        0.059605849715634,  0.329505331313828,  0.501692882467059,
        0.000000000000000, -0.218104134283380,  0.000000000000000,
        0.000000000000000, -0.391123811866453,  0.000000000000000,
        0.000000000000000,  0.000000000000000, -0.062544706021092,
        0.000000000000000, -0.108330608572991,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.501520641074811,  0.000000000000000,  0.000000000000000,
        1.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.147927358478970,  0.309505179870655,
       -0.100988447811255, -0.083880016685293,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.089333516637993,  0.000000000000000,
        0.000000000000000,  0.342899055855409,  0.000000000000000,
       -0.180410116244944,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.501520641074811,  0.000000000000000,  0.000000000000000,
        1.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.078715737682824,  0.164695217981853,  0.529105446483568,
        0.439469807103769,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.062253543338693,  0.000000000000000,
        0.000000000000000,  0.180410116244944,  0.000000000000000,
        0.000000000000000, -0.224124179944020,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        1.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.010755059186520,
        0.021770630621996,  0.043677731061311,  0.008151460212183,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.062253543338693,
        0.000000000000000,  0.000000000000000,  0.180410116244944,
        0.000000000000000, -0.224124179944020,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        1.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.025963452940715,  0.052555800376472,
       -0.105441048282554, -0.019678185860883, -0.006457144066568,
       -0.070665963942989, -0.076197708983921,  0.000000000000000,
        0.095489194679018,  0.000000000000000,  0.000000000000000,
        0.062544706021092,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.093891255171399,  0.000000000000000,
        0.089962792984985,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        1.000000000000001,  0.000000000000000,  0.000000000000000,
       -0.006792877136544, -0.013750293367576,  0.285773499802488,
        0.053333157577861,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.041951212424467,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        1.000000000000000,  0.000000000000000,  0.005723028021282,
        0.011584680932934, -0.228839296924101, -0.042707676853505,
       -0.011184101595088, -0.122397039915088, -0.131978303380498,
        0.000000000000000,  0.165392136757895,  0.000000000000000,
        0.000000000000000,  0.108330608572991,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.089962792984985,
        0.000000000000000,  0.197771340665263,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        0.000000000000000,  0.000000000000000,  0.000000000000000,
        1.000000000000000, -0.023210693813002, -0.046983603969048,
        0.037334435775198,  0.006967627675958,  0.033253425203908,
        0.240404328848086,  0.420483682453648,  0.105654271868166,
        0.048681748856291,  0.255056682446983,  0.219175793213426,
        0.100988447811255,  0.529105446483568,  0.043677731061311,
        0.105441048282554,  0.285773499802488,  0.228839296924101,
        0.037334435775198,  0.000104734953463,  0.005168519068412,
        0.063349060947094,  0.001871814694006,  0.008491793599620,
        0.004518689472437,  0.032607081177624,  0.147927358478970,
        0.078715737682824,  0.010755059186520,  0.025963452940715,
       -0.006792877136544,  0.005723028021282, -0.023210693813002,
        1.000000000000000,  0.658292049339330,  0.008118968797873,
        0.107356188553075,  0.068076201809780,  0.375934821334508,
        0.678093402858836,  0.055743215005449,  0.025684500449954,
        0.134567956759564,  0.182045269436310,  0.083880016685293,
        0.439469807103769,  0.008151460212183,  0.019678185860883,
        0.053333157577861,  0.042707676853505,  0.006967627675958,
        0.014626666566694,  0.098291977680733,  0.254121532334353,
        0.014195841148136,  0.064401755894416,  0.034269737359089,
        0.068223083469529,  0.309505179870655,  0.164695217981853,
        0.021770630621996,  0.052555800376472, -0.013750293367576,
        0.011584680932934, -0.046983603969048,  0.658292049339330,
        1.000000000000001,  0.107356188553075,  0.327910319525025,
        0.000104734953463,  0.005168519068412,  0.063349060947094,
       -0.001871814694006, -0.008491793599620,  0.004518689472437,
       -0.032607081177624, -0.147927358478970,  0.078715737682824,
        0.010755059186520, -0.025963452940715, -0.006792877136544,
       -0.005723028021282, -0.023210693813002,  0.033253425203908,
        0.240404328848086,  0.420483682453648, -0.105654271868166,
       -0.048681748856291,  0.255056682446983, -0.219175793213426,
       -0.100988447811255,  0.529105446483568,  0.043677731061311,
       -0.105441048282554,  0.285773499802488, -0.228839296924101,
        0.037334435775198,  0.008118968797873,  0.107356188553075,
        1.000000000000000,  0.658292049339330,  0.014626666566694,
        0.098291977680733,  0.254121532334353, -0.014195841148136,
       -0.064401755894416,  0.034269737359089, -0.068223083469529,
       -0.309505179870655,  0.164695217981853,  0.021770630621996,
       -0.052555800376472, -0.013750293367576, -0.011584680932934,
       -0.046983603969048,  0.068076201809780,  0.375934821334508,
        0.678093402858836, -0.055743215005449, -0.025684500449954,
        0.134567956759564, -0.182045269436310, -0.083880016685293,
        0.439469807103769,  0.008151460212183, -0.019678185860883,
        0.053333157577861, -0.042707676853505,  0.006967627675958,
        0.107356188553075,  0.327910319525025,  0.658292049339330,
        1.000000000000001
};