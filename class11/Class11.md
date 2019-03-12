Class 11: Structural Bioinformatics 1
================

Section 1: PDB
--------------

Q1. Determine the percentage of structures solved by X-Ray and Electron Microscopy. What proportion of structures are protein?

``` r
stats <- read.csv("Data Export Summary.csv", row.names = 1)
stats
```

    ##                     Proteins Nucleic.Acids Protein.NA.Complex Other  Total
    ## X-Ray                 124770          1993               6451    10 133224
    ## NMR                    10988          1273                257     8  12526
    ## Electron Microscopy     2057            31                723     0   2811
    ## Other                    250             4                  6    13    273
    ## Multi Method             127             5                  2     1    135

``` r
total <- sum(stats$Total)
stats$Total[1]/total *100
```

    ## [1] 89.43069

``` r
stats$Total[3]/total *100
```

    ## [1] 1.88697

``` r
percent.by.method <- stats$Total/total*100
names(percent.by.method) <rownames(stats)
```

    ## logical(0)

``` r
percent.by.method
```

    ## [1] 89.43068692  8.40846082  1.88696977  0.18325960  0.09062288

``` r
sum(stats$Proteins)/total*100
```

    ## [1] 92.76561

Section 3

``` r
library(bio3d)

pdb <- read.pdb("1hsg")
```

    ##   Note: Accessing on-line PDB file

``` r
pdb$atom$resid
```

    ##    [1] "PRO" "PRO" "PRO" "PRO" "PRO" "PRO" "PRO" "GLN" "GLN" "GLN" "GLN"
    ##   [12] "GLN" "GLN" "GLN" "GLN" "GLN" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE"
    ##   [23] "ILE" "ILE" "THR" "THR" "THR" "THR" "THR" "THR" "THR" "LEU" "LEU"
    ##   [34] "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "TRP" "TRP" "TRP" "TRP" "TRP"
    ##   [45] "TRP" "TRP" "TRP" "TRP" "TRP" "TRP" "TRP" "TRP" "TRP" "GLN" "GLN"
    ##   [56] "GLN" "GLN" "GLN" "GLN" "GLN" "GLN" "GLN" "ARG" "ARG" "ARG" "ARG"
    ##   [67] "ARG" "ARG" "ARG" "ARG" "ARG" "ARG" "ARG" "PRO" "PRO" "PRO" "PRO"
    ##   [78] "PRO" "PRO" "PRO" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU"
    ##   [89] "VAL" "VAL" "VAL" "VAL" "VAL" "VAL" "VAL" "THR" "THR" "THR" "THR"
    ##  [100] "THR" "THR" "THR" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE"
    ##  [111] "LYS" "LYS" "LYS" "LYS" "LYS" "LYS" "LYS" "LYS" "LYS" "ILE" "ILE"
    ##  [122] "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "GLY" "GLY" "GLY" "GLY" "GLY"
    ##  [133] "GLY" "GLY" "GLY" "GLN" "GLN" "GLN" "GLN" "GLN" "GLN" "GLN" "GLN"
    ##  [144] "GLN" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LYS" "LYS"
    ##  [155] "LYS" "LYS" "LYS" "LYS" "LYS" "LYS" "LYS" "GLU" "GLU" "GLU" "GLU"
    ##  [166] "GLU" "GLU" "GLU" "GLU" "GLU" "ALA" "ALA" "ALA" "ALA" "ALA" "LEU"
    ##  [177] "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU"
    ##  [188] "LEU" "LEU" "LEU" "LEU" "ASP" "ASP" "ASP" "ASP" "ASP" "ASP" "ASP"
    ##  [199] "ASP" "THR" "THR" "THR" "THR" "THR" "THR" "THR" "GLY" "GLY" "GLY"
    ##  [210] "GLY" "ALA" "ALA" "ALA" "ALA" "ALA" "ASP" "ASP" "ASP" "ASP" "ASP"
    ##  [221] "ASP" "ASP" "ASP" "ASP" "ASP" "ASP" "ASP" "ASP" "ASP" "ASP" "ASP"
    ##  [232] "THR" "THR" "THR" "THR" "THR" "THR" "THR" "VAL" "VAL" "VAL" "VAL"
    ##  [243] "VAL" "VAL" "VAL" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU"
    ##  [254] "GLU" "GLU" "GLU" "GLU" "GLU" "GLU" "GLU" "GLU" "GLU" "GLU" "GLU"
    ##  [265] "GLU" "GLU" "GLU" "GLU" "GLU" "GLU" "GLU" "MET" "MET" "MET" "MET"
    ##  [276] "MET" "MET" "MET" "MET" "SER" "SER" "SER" "SER" "SER" "SER" "LEU"
    ##  [287] "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "PRO" "PRO" "PRO" "PRO"
    ##  [298] "PRO" "PRO" "PRO" "GLY" "GLY" "GLY" "GLY" "ARG" "ARG" "ARG" "ARG"
    ##  [309] "ARG" "ARG" "ARG" "ARG" "ARG" "ARG" "ARG" "TRP" "TRP" "TRP" "TRP"
    ##  [320] "TRP" "TRP" "TRP" "TRP" "TRP" "TRP" "TRP" "TRP" "TRP" "TRP" "LYS"
    ##  [331] "LYS" "LYS" "LYS" "LYS" "LYS" "LYS" "LYS" "LYS" "PRO" "PRO" "PRO"
    ##  [342] "PRO" "PRO" "PRO" "PRO" "LYS" "LYS" "LYS" "LYS" "LYS" "LYS" "LYS"
    ##  [353] "LYS" "LYS" "MET" "MET" "MET" "MET" "MET" "MET" "MET" "MET" "ILE"
    ##  [364] "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "GLY" "GLY" "GLY" "GLY"
    ##  [375] "GLY" "GLY" "GLY" "GLY" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE"
    ##  [386] "ILE" "GLY" "GLY" "GLY" "GLY" "GLY" "GLY" "GLY" "GLY" "PHE" "PHE"
    ##  [397] "PHE" "PHE" "PHE" "PHE" "PHE" "PHE" "PHE" "PHE" "PHE" "ILE" "ILE"
    ##  [408] "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "LYS" "LYS" "LYS" "LYS" "LYS"
    ##  [419] "LYS" "LYS" "LYS" "LYS" "VAL" "VAL" "VAL" "VAL" "VAL" "VAL" "VAL"
    ##  [430] "ARG" "ARG" "ARG" "ARG" "ARG" "ARG" "ARG" "ARG" "ARG" "ARG" "ARG"
    ##  [441] "GLN" "GLN" "GLN" "GLN" "GLN" "GLN" "GLN" "GLN" "GLN" "TYR" "TYR"
    ##  [452] "TYR" "TYR" "TYR" "TYR" "TYR" "TYR" "TYR" "TYR" "TYR" "TYR" "ASP"
    ##  [463] "ASP" "ASP" "ASP" "ASP" "ASP" "ASP" "ASP" "GLN" "GLN" "GLN" "GLN"
    ##  [474] "GLN" "GLN" "GLN" "GLN" "GLN" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE"
    ##  [485] "ILE" "ILE" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "ILE"
    ##  [496] "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "GLU" "GLU" "GLU" "GLU"
    ##  [507] "GLU" "GLU" "GLU" "GLU" "GLU" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE"
    ##  [518] "ILE" "ILE" "CYS" "CYS" "CYS" "CYS" "CYS" "CYS" "GLY" "GLY" "GLY"
    ##  [529] "GLY" "HIS" "HIS" "HIS" "HIS" "HIS" "HIS" "HIS" "HIS" "HIS" "HIS"
    ##  [540] "LYS" "LYS" "LYS" "LYS" "LYS" "LYS" "LYS" "LYS" "LYS" "ALA" "ALA"
    ##  [551] "ALA" "ALA" "ALA" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE"
    ##  [562] "GLY" "GLY" "GLY" "GLY" "THR" "THR" "THR" "THR" "THR" "THR" "THR"
    ##  [573] "VAL" "VAL" "VAL" "VAL" "VAL" "VAL" "VAL" "LEU" "LEU" "LEU" "LEU"
    ##  [584] "LEU" "LEU" "LEU" "LEU" "VAL" "VAL" "VAL" "VAL" "VAL" "VAL" "VAL"
    ##  [595] "GLY" "GLY" "GLY" "GLY" "PRO" "PRO" "PRO" "PRO" "PRO" "PRO" "PRO"
    ##  [606] "THR" "THR" "THR" "THR" "THR" "THR" "THR" "PRO" "PRO" "PRO" "PRO"
    ##  [617] "PRO" "PRO" "PRO" "VAL" "VAL" "VAL" "VAL" "VAL" "VAL" "VAL" "ASN"
    ##  [628] "ASN" "ASN" "ASN" "ASN" "ASN" "ASN" "ASN" "ILE" "ILE" "ILE" "ILE"
    ##  [639] "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE"
    ##  [650] "ILE" "GLY" "GLY" "GLY" "GLY" "ARG" "ARG" "ARG" "ARG" "ARG" "ARG"
    ##  [661] "ARG" "ARG" "ARG" "ARG" "ARG" "ASN" "ASN" "ASN" "ASN" "ASN" "ASN"
    ##  [672] "ASN" "ASN" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU"
    ##  [683] "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "THR" "THR" "THR" "THR"
    ##  [694] "THR" "THR" "THR" "GLN" "GLN" "GLN" "GLN" "GLN" "GLN" "GLN" "GLN"
    ##  [705] "GLN" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "GLY" "GLY"
    ##  [716] "GLY" "GLY" "CYS" "CYS" "CYS" "CYS" "CYS" "CYS" "THR" "THR" "THR"
    ##  [727] "THR" "THR" "THR" "THR" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU"
    ##  [738] "LEU" "ASN" "ASN" "ASN" "ASN" "ASN" "ASN" "ASN" "ASN" "PHE" "PHE"
    ##  [749] "PHE" "PHE" "PHE" "PHE" "PHE" "PHE" "PHE" "PHE" "PHE" "PRO" "PRO"
    ##  [760] "PRO" "PRO" "PRO" "PRO" "PRO" "GLN" "GLN" "GLN" "GLN" "GLN" "GLN"
    ##  [771] "GLN" "GLN" "GLN" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE"
    ##  [782] "THR" "THR" "THR" "THR" "THR" "THR" "THR" "LEU" "LEU" "LEU" "LEU"
    ##  [793] "LEU" "LEU" "LEU" "LEU" "TRP" "TRP" "TRP" "TRP" "TRP" "TRP" "TRP"
    ##  [804] "TRP" "TRP" "TRP" "TRP" "TRP" "TRP" "TRP" "GLN" "GLN" "GLN" "GLN"
    ##  [815] "GLN" "GLN" "GLN" "GLN" "GLN" "ARG" "ARG" "ARG" "ARG" "ARG" "ARG"
    ##  [826] "ARG" "ARG" "ARG" "ARG" "ARG" "PRO" "PRO" "PRO" "PRO" "PRO" "PRO"
    ##  [837] "PRO" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "VAL" "VAL"
    ##  [848] "VAL" "VAL" "VAL" "VAL" "VAL" "THR" "THR" "THR" "THR" "THR" "THR"
    ##  [859] "THR" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "LYS" "LYS"
    ##  [870] "LYS" "LYS" "LYS" "LYS" "LYS" "LYS" "LYS" "ILE" "ILE" "ILE" "ILE"
    ##  [881] "ILE" "ILE" "ILE" "ILE" "GLY" "GLY" "GLY" "GLY" "GLY" "GLY" "GLY"
    ##  [892] "GLY" "GLN" "GLN" "GLN" "GLN" "GLN" "GLN" "GLN" "GLN" "GLN" "LEU"
    ##  [903] "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LYS" "LYS" "LYS" "LYS"
    ##  [914] "LYS" "LYS" "LYS" "LYS" "LYS" "GLU" "GLU" "GLU" "GLU" "GLU" "GLU"
    ##  [925] "GLU" "GLU" "GLU" "ALA" "ALA" "ALA" "ALA" "ALA" "LEU" "LEU" "LEU"
    ##  [936] "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU"
    ##  [947] "LEU" "LEU" "ASP" "ASP" "ASP" "ASP" "ASP" "ASP" "ASP" "ASP" "THR"
    ##  [958] "THR" "THR" "THR" "THR" "THR" "THR" "GLY" "GLY" "GLY" "GLY" "ALA"
    ##  [969] "ALA" "ALA" "ALA" "ALA" "ASP" "ASP" "ASP" "ASP" "ASP" "ASP" "ASP"
    ##  [980] "ASP" "ASP" "ASP" "ASP" "ASP" "ASP" "ASP" "ASP" "ASP" "THR" "THR"
    ##  [991] "THR" "THR" "THR" "THR" "THR" "VAL" "VAL" "VAL" "VAL" "VAL" "VAL"
    ## [1002] "VAL" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "GLU" "GLU"
    ## [1013] "GLU" "GLU" "GLU" "GLU" "GLU" "GLU" "GLU" "GLU" "GLU" "GLU" "GLU"
    ## [1024] "GLU" "GLU" "GLU" "GLU" "GLU" "MET" "MET" "MET" "MET" "MET" "MET"
    ## [1035] "MET" "MET" "SER" "SER" "SER" "SER" "SER" "SER" "LEU" "LEU" "LEU"
    ## [1046] "LEU" "LEU" "LEU" "LEU" "LEU" "PRO" "PRO" "PRO" "PRO" "PRO" "PRO"
    ## [1057] "PRO" "GLY" "GLY" "GLY" "GLY" "ARG" "ARG" "ARG" "ARG" "ARG" "ARG"
    ## [1068] "ARG" "ARG" "ARG" "ARG" "ARG" "TRP" "TRP" "TRP" "TRP" "TRP" "TRP"
    ## [1079] "TRP" "TRP" "TRP" "TRP" "TRP" "TRP" "TRP" "TRP" "LYS" "LYS" "LYS"
    ## [1090] "LYS" "LYS" "LYS" "LYS" "LYS" "LYS" "PRO" "PRO" "PRO" "PRO" "PRO"
    ## [1101] "PRO" "PRO" "LYS" "LYS" "LYS" "LYS" "LYS" "LYS" "LYS" "LYS" "LYS"
    ## [1112] "MET" "MET" "MET" "MET" "MET" "MET" "MET" "MET" "ILE" "ILE" "ILE"
    ## [1123] "ILE" "ILE" "ILE" "ILE" "ILE" "GLY" "GLY" "GLY" "GLY" "GLY" "GLY"
    ## [1134] "GLY" "GLY" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "GLY"
    ## [1145] "GLY" "GLY" "GLY" "GLY" "GLY" "GLY" "GLY" "PHE" "PHE" "PHE" "PHE"
    ## [1156] "PHE" "PHE" "PHE" "PHE" "PHE" "PHE" "PHE" "ILE" "ILE" "ILE" "ILE"
    ## [1167] "ILE" "ILE" "ILE" "ILE" "LYS" "LYS" "LYS" "LYS" "LYS" "LYS" "LYS"
    ## [1178] "LYS" "LYS" "VAL" "VAL" "VAL" "VAL" "VAL" "VAL" "VAL" "ARG" "ARG"
    ## [1189] "ARG" "ARG" "ARG" "ARG" "ARG" "ARG" "ARG" "ARG" "ARG" "GLN" "GLN"
    ## [1200] "GLN" "GLN" "GLN" "GLN" "GLN" "GLN" "GLN" "TYR" "TYR" "TYR" "TYR"
    ## [1211] "TYR" "TYR" "TYR" "TYR" "TYR" "TYR" "TYR" "TYR" "ASP" "ASP" "ASP"
    ## [1222] "ASP" "ASP" "ASP" "ASP" "ASP" "GLN" "GLN" "GLN" "GLN" "GLN" "GLN"
    ## [1233] "GLN" "GLN" "GLN" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE"
    ## [1244] "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "ILE" "ILE" "ILE"
    ## [1255] "ILE" "ILE" "ILE" "ILE" "ILE" "GLU" "GLU" "GLU" "GLU" "GLU" "GLU"
    ## [1266] "GLU" "GLU" "GLU" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE"
    ## [1277] "CYS" "CYS" "CYS" "CYS" "CYS" "CYS" "GLY" "GLY" "GLY" "GLY" "HIS"
    ## [1288] "HIS" "HIS" "HIS" "HIS" "HIS" "HIS" "HIS" "HIS" "HIS" "LYS" "LYS"
    ## [1299] "LYS" "LYS" "LYS" "LYS" "LYS" "LYS" "LYS" "ALA" "ALA" "ALA" "ALA"
    ## [1310] "ALA" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "GLY" "GLY"
    ## [1321] "GLY" "GLY" "THR" "THR" "THR" "THR" "THR" "THR" "THR" "VAL" "VAL"
    ## [1332] "VAL" "VAL" "VAL" "VAL" "VAL" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU"
    ## [1343] "LEU" "LEU" "VAL" "VAL" "VAL" "VAL" "VAL" "VAL" "VAL" "GLY" "GLY"
    ## [1354] "GLY" "GLY" "PRO" "PRO" "PRO" "PRO" "PRO" "PRO" "PRO" "THR" "THR"
    ## [1365] "THR" "THR" "THR" "THR" "THR" "PRO" "PRO" "PRO" "PRO" "PRO" "PRO"
    ## [1376] "PRO" "VAL" "VAL" "VAL" "VAL" "VAL" "VAL" "VAL" "ASN" "ASN" "ASN"
    ## [1387] "ASN" "ASN" "ASN" "ASN" "ASN" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE"
    ## [1398] "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "GLY"
    ## [1409] "GLY" "GLY" "GLY" "ARG" "ARG" "ARG" "ARG" "ARG" "ARG" "ARG" "ARG"
    ## [1420] "ARG" "ARG" "ARG" "ASN" "ASN" "ASN" "ASN" "ASN" "ASN" "ASN" "ASN"
    ## [1431] "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU"
    ## [1442] "LEU" "LEU" "LEU" "LEU" "LEU" "THR" "THR" "THR" "THR" "THR" "THR"
    ## [1453] "THR" "GLN" "GLN" "GLN" "GLN" "GLN" "GLN" "GLN" "GLN" "GLN" "ILE"
    ## [1464] "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "GLY" "GLY" "GLY" "GLY"
    ## [1475] "CYS" "CYS" "CYS" "CYS" "CYS" "CYS" "THR" "THR" "THR" "THR" "THR"
    ## [1486] "THR" "THR" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "ASN"
    ## [1497] "ASN" "ASN" "ASN" "ASN" "ASN" "ASN" "ASN" "PHE" "PHE" "PHE" "PHE"
    ## [1508] "PHE" "PHE" "PHE" "PHE" "PHE" "PHE" "PHE" "MK1" "MK1" "MK1" "MK1"
    ## [1519] "MK1" "MK1" "MK1" "MK1" "MK1" "MK1" "MK1" "MK1" "MK1" "MK1" "MK1"
    ## [1530] "MK1" "MK1" "MK1" "MK1" "MK1" "MK1" "MK1" "MK1" "MK1" "MK1" "MK1"
    ## [1541] "MK1" "MK1" "MK1" "MK1" "MK1" "MK1" "MK1" "MK1" "MK1" "MK1" "MK1"
    ## [1552] "MK1" "MK1" "MK1" "MK1" "MK1" "MK1" "MK1" "MK1" "HOH" "HOH" "HOH"
    ## [1563] "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH"
    ## [1574] "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH"
    ## [1585] "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH"
    ## [1596] "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH"
    ## [1607] "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH"
    ## [1618] "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH"
    ## [1629] "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH"
    ## [1640] "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH"
    ## [1651] "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH"
    ## [1662] "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH"
    ## [1673] "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH"
    ## [1684] "HOH" "HOH" "HOH"

``` r
aa321(pdb$atom$resid)
```

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ##    [1] "P" "P" "P" "P" "P" "P" "P" "Q" "Q" "Q" "Q" "Q" "Q" "Q" "Q" "Q" "I"
    ##   [18] "I" "I" "I" "I" "I" "I" "I" "T" "T" "T" "T" "T" "T" "T" "L" "L" "L"
    ##   [35] "L" "L" "L" "L" "L" "W" "W" "W" "W" "W" "W" "W" "W" "W" "W" "W" "W"
    ##   [52] "W" "W" "Q" "Q" "Q" "Q" "Q" "Q" "Q" "Q" "Q" "R" "R" "R" "R" "R" "R"
    ##   [69] "R" "R" "R" "R" "R" "P" "P" "P" "P" "P" "P" "P" "L" "L" "L" "L" "L"
    ##   [86] "L" "L" "L" "V" "V" "V" "V" "V" "V" "V" "T" "T" "T" "T" "T" "T" "T"
    ##  [103] "I" "I" "I" "I" "I" "I" "I" "I" "K" "K" "K" "K" "K" "K" "K" "K" "K"
    ##  [120] "I" "I" "I" "I" "I" "I" "I" "I" "G" "G" "G" "G" "G" "G" "G" "G" "Q"
    ##  [137] "Q" "Q" "Q" "Q" "Q" "Q" "Q" "Q" "L" "L" "L" "L" "L" "L" "L" "L" "K"
    ##  [154] "K" "K" "K" "K" "K" "K" "K" "K" "E" "E" "E" "E" "E" "E" "E" "E" "E"
    ##  [171] "A" "A" "A" "A" "A" "L" "L" "L" "L" "L" "L" "L" "L" "L" "L" "L" "L"
    ##  [188] "L" "L" "L" "L" "D" "D" "D" "D" "D" "D" "D" "D" "T" "T" "T" "T" "T"
    ##  [205] "T" "T" "G" "G" "G" "G" "A" "A" "A" "A" "A" "D" "D" "D" "D" "D" "D"
    ##  [222] "D" "D" "D" "D" "D" "D" "D" "D" "D" "D" "T" "T" "T" "T" "T" "T" "T"
    ##  [239] "V" "V" "V" "V" "V" "V" "V" "L" "L" "L" "L" "L" "L" "L" "L" "E" "E"
    ##  [256] "E" "E" "E" "E" "E" "E" "E" "E" "E" "E" "E" "E" "E" "E" "E" "E" "M"
    ##  [273] "M" "M" "M" "M" "M" "M" "M" "S" "S" "S" "S" "S" "S" "L" "L" "L" "L"
    ##  [290] "L" "L" "L" "L" "P" "P" "P" "P" "P" "P" "P" "G" "G" "G" "G" "R" "R"
    ##  [307] "R" "R" "R" "R" "R" "R" "R" "R" "R" "W" "W" "W" "W" "W" "W" "W" "W"
    ##  [324] "W" "W" "W" "W" "W" "W" "K" "K" "K" "K" "K" "K" "K" "K" "K" "P" "P"
    ##  [341] "P" "P" "P" "P" "P" "K" "K" "K" "K" "K" "K" "K" "K" "K" "M" "M" "M"
    ##  [358] "M" "M" "M" "M" "M" "I" "I" "I" "I" "I" "I" "I" "I" "G" "G" "G" "G"
    ##  [375] "G" "G" "G" "G" "I" "I" "I" "I" "I" "I" "I" "I" "G" "G" "G" "G" "G"
    ##  [392] "G" "G" "G" "F" "F" "F" "F" "F" "F" "F" "F" "F" "F" "F" "I" "I" "I"
    ##  [409] "I" "I" "I" "I" "I" "K" "K" "K" "K" "K" "K" "K" "K" "K" "V" "V" "V"
    ##  [426] "V" "V" "V" "V" "R" "R" "R" "R" "R" "R" "R" "R" "R" "R" "R" "Q" "Q"
    ##  [443] "Q" "Q" "Q" "Q" "Q" "Q" "Q" "Y" "Y" "Y" "Y" "Y" "Y" "Y" "Y" "Y" "Y"
    ##  [460] "Y" "Y" "D" "D" "D" "D" "D" "D" "D" "D" "Q" "Q" "Q" "Q" "Q" "Q" "Q"
    ##  [477] "Q" "Q" "I" "I" "I" "I" "I" "I" "I" "I" "L" "L" "L" "L" "L" "L" "L"
    ##  [494] "L" "I" "I" "I" "I" "I" "I" "I" "I" "E" "E" "E" "E" "E" "E" "E" "E"
    ##  [511] "E" "I" "I" "I" "I" "I" "I" "I" "I" "C" "C" "C" "C" "C" "C" "G" "G"
    ##  [528] "G" "G" "H" "H" "H" "H" "H" "H" "H" "H" "H" "H" "K" "K" "K" "K" "K"
    ##  [545] "K" "K" "K" "K" "A" "A" "A" "A" "A" "I" "I" "I" "I" "I" "I" "I" "I"
    ##  [562] "G" "G" "G" "G" "T" "T" "T" "T" "T" "T" "T" "V" "V" "V" "V" "V" "V"
    ##  [579] "V" "L" "L" "L" "L" "L" "L" "L" "L" "V" "V" "V" "V" "V" "V" "V" "G"
    ##  [596] "G" "G" "G" "P" "P" "P" "P" "P" "P" "P" "T" "T" "T" "T" "T" "T" "T"
    ##  [613] "P" "P" "P" "P" "P" "P" "P" "V" "V" "V" "V" "V" "V" "V" "N" "N" "N"
    ##  [630] "N" "N" "N" "N" "N" "I" "I" "I" "I" "I" "I" "I" "I" "I" "I" "I" "I"
    ##  [647] "I" "I" "I" "I" "G" "G" "G" "G" "R" "R" "R" "R" "R" "R" "R" "R" "R"
    ##  [664] "R" "R" "N" "N" "N" "N" "N" "N" "N" "N" "L" "L" "L" "L" "L" "L" "L"
    ##  [681] "L" "L" "L" "L" "L" "L" "L" "L" "L" "T" "T" "T" "T" "T" "T" "T" "Q"
    ##  [698] "Q" "Q" "Q" "Q" "Q" "Q" "Q" "Q" "I" "I" "I" "I" "I" "I" "I" "I" "G"
    ##  [715] "G" "G" "G" "C" "C" "C" "C" "C" "C" "T" "T" "T" "T" "T" "T" "T" "L"
    ##  [732] "L" "L" "L" "L" "L" "L" "L" "N" "N" "N" "N" "N" "N" "N" "N" "F" "F"
    ##  [749] "F" "F" "F" "F" "F" "F" "F" "F" "F" "P" "P" "P" "P" "P" "P" "P" "Q"
    ##  [766] "Q" "Q" "Q" "Q" "Q" "Q" "Q" "Q" "I" "I" "I" "I" "I" "I" "I" "I" "T"
    ##  [783] "T" "T" "T" "T" "T" "T" "L" "L" "L" "L" "L" "L" "L" "L" "W" "W" "W"
    ##  [800] "W" "W" "W" "W" "W" "W" "W" "W" "W" "W" "W" "Q" "Q" "Q" "Q" "Q" "Q"
    ##  [817] "Q" "Q" "Q" "R" "R" "R" "R" "R" "R" "R" "R" "R" "R" "R" "P" "P" "P"
    ##  [834] "P" "P" "P" "P" "L" "L" "L" "L" "L" "L" "L" "L" "V" "V" "V" "V" "V"
    ##  [851] "V" "V" "T" "T" "T" "T" "T" "T" "T" "I" "I" "I" "I" "I" "I" "I" "I"
    ##  [868] "K" "K" "K" "K" "K" "K" "K" "K" "K" "I" "I" "I" "I" "I" "I" "I" "I"
    ##  [885] "G" "G" "G" "G" "G" "G" "G" "G" "Q" "Q" "Q" "Q" "Q" "Q" "Q" "Q" "Q"
    ##  [902] "L" "L" "L" "L" "L" "L" "L" "L" "K" "K" "K" "K" "K" "K" "K" "K" "K"
    ##  [919] "E" "E" "E" "E" "E" "E" "E" "E" "E" "A" "A" "A" "A" "A" "L" "L" "L"
    ##  [936] "L" "L" "L" "L" "L" "L" "L" "L" "L" "L" "L" "L" "L" "D" "D" "D" "D"
    ##  [953] "D" "D" "D" "D" "T" "T" "T" "T" "T" "T" "T" "G" "G" "G" "G" "A" "A"
    ##  [970] "A" "A" "A" "D" "D" "D" "D" "D" "D" "D" "D" "D" "D" "D" "D" "D" "D"
    ##  [987] "D" "D" "T" "T" "T" "T" "T" "T" "T" "V" "V" "V" "V" "V" "V" "V" "L"
    ## [1004] "L" "L" "L" "L" "L" "L" "L" "E" "E" "E" "E" "E" "E" "E" "E" "E" "E"
    ## [1021] "E" "E" "E" "E" "E" "E" "E" "E" "M" "M" "M" "M" "M" "M" "M" "M" "S"
    ## [1038] "S" "S" "S" "S" "S" "L" "L" "L" "L" "L" "L" "L" "L" "P" "P" "P" "P"
    ## [1055] "P" "P" "P" "G" "G" "G" "G" "R" "R" "R" "R" "R" "R" "R" "R" "R" "R"
    ## [1072] "R" "W" "W" "W" "W" "W" "W" "W" "W" "W" "W" "W" "W" "W" "W" "K" "K"
    ## [1089] "K" "K" "K" "K" "K" "K" "K" "P" "P" "P" "P" "P" "P" "P" "K" "K" "K"
    ## [1106] "K" "K" "K" "K" "K" "K" "M" "M" "M" "M" "M" "M" "M" "M" "I" "I" "I"
    ## [1123] "I" "I" "I" "I" "I" "G" "G" "G" "G" "G" "G" "G" "G" "I" "I" "I" "I"
    ## [1140] "I" "I" "I" "I" "G" "G" "G" "G" "G" "G" "G" "G" "F" "F" "F" "F" "F"
    ## [1157] "F" "F" "F" "F" "F" "F" "I" "I" "I" "I" "I" "I" "I" "I" "K" "K" "K"
    ## [1174] "K" "K" "K" "K" "K" "K" "V" "V" "V" "V" "V" "V" "V" "R" "R" "R" "R"
    ## [1191] "R" "R" "R" "R" "R" "R" "R" "Q" "Q" "Q" "Q" "Q" "Q" "Q" "Q" "Q" "Y"
    ## [1208] "Y" "Y" "Y" "Y" "Y" "Y" "Y" "Y" "Y" "Y" "Y" "D" "D" "D" "D" "D" "D"
    ## [1225] "D" "D" "Q" "Q" "Q" "Q" "Q" "Q" "Q" "Q" "Q" "I" "I" "I" "I" "I" "I"
    ## [1242] "I" "I" "L" "L" "L" "L" "L" "L" "L" "L" "I" "I" "I" "I" "I" "I" "I"
    ## [1259] "I" "E" "E" "E" "E" "E" "E" "E" "E" "E" "I" "I" "I" "I" "I" "I" "I"
    ## [1276] "I" "C" "C" "C" "C" "C" "C" "G" "G" "G" "G" "H" "H" "H" "H" "H" "H"
    ## [1293] "H" "H" "H" "H" "K" "K" "K" "K" "K" "K" "K" "K" "K" "A" "A" "A" "A"
    ## [1310] "A" "I" "I" "I" "I" "I" "I" "I" "I" "G" "G" "G" "G" "T" "T" "T" "T"
    ## [1327] "T" "T" "T" "V" "V" "V" "V" "V" "V" "V" "L" "L" "L" "L" "L" "L" "L"
    ## [1344] "L" "V" "V" "V" "V" "V" "V" "V" "G" "G" "G" "G" "P" "P" "P" "P" "P"
    ## [1361] "P" "P" "T" "T" "T" "T" "T" "T" "T" "P" "P" "P" "P" "P" "P" "P" "V"
    ## [1378] "V" "V" "V" "V" "V" "V" "N" "N" "N" "N" "N" "N" "N" "N" "I" "I" "I"
    ## [1395] "I" "I" "I" "I" "I" "I" "I" "I" "I" "I" "I" "I" "I" "G" "G" "G" "G"
    ## [1412] "R" "R" "R" "R" "R" "R" "R" "R" "R" "R" "R" "N" "N" "N" "N" "N" "N"
    ## [1429] "N" "N" "L" "L" "L" "L" "L" "L" "L" "L" "L" "L" "L" "L" "L" "L" "L"
    ## [1446] "L" "T" "T" "T" "T" "T" "T" "T" "Q" "Q" "Q" "Q" "Q" "Q" "Q" "Q" "Q"
    ## [1463] "I" "I" "I" "I" "I" "I" "I" "I" "G" "G" "G" "G" "C" "C" "C" "C" "C"
    ## [1480] "C" "T" "T" "T" "T" "T" "T" "T" "L" "L" "L" "L" "L" "L" "L" "L" "N"
    ## [1497] "N" "N" "N" "N" "N" "N" "N" "F" "F" "F" "F" "F" "F" "F" "F" "F" "F"
    ## [1514] "F" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X"
    ## [1531] "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X"
    ## [1548] "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X"
    ## [1565] "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X"
    ## [1582] "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X"
    ## [1599] "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X"
    ## [1616] "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X"
    ## [1633] "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X"
    ## [1650] "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X"
    ## [1667] "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X"
    ## [1684] "X" "X" "X"

Q8. Use the Bio3D write.pdb() function to write out a protein only PDB file for viewing in VMD. We want to select out the protein and drug only parts of these molecular PDB files.

``` r
pdb <- read.pdb("1hsg")
```

    ##   Note: Accessing on-line PDB file

    ## Warning in get.pdb(file, path = tempdir(), verbose = FALSE): /var/folders/
    ## 60/y7f0bdpn5jb9zcc5496qgb9h0000gn/T//RtmpCLSMT1/1hsg.pdb exists. Skipping
    ## download

``` r
prot.inds <- atom.select(pdb, "protein")
prot.inds
```

    ## 
    ##  Call:  atom.select.pdb(pdb = pdb, string = "protein")
    ## 
    ##    Atom Indices#: 1514  ($atom)
    ##    XYZ  Indices#: 4542  ($xyz)
    ## 
    ## + attr: atom, xyz, call

``` r
prot.pdb <- trim.pdb(pdb, prot.inds)
write.pdb(prot.pdb, file = "protein.pdb")
```

``` r
lig.inds <- atom.select(pdb, "ligand")
lig.pdb <- trim.pdb(pdb, lig.inds)
write.pdb(prot.pdb, file = "ligand.pdb")
```
