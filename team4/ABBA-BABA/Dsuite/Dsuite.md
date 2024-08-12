```sh
Usage:
a) Dsuite Dtrios SNP.bg.vcf.gz  SETS.txt
b) Dsuite Dquartets SNP.bg.vcf.gz SETS.txt
```

 Population/species map SETS.txt: a text file with one individual per row and a tab separating the individual’s name from the name of the species/population it belongs to, as shown below:

```
A01     Malus_domestica2
A02     Malus_domestica2
A03     Malus_domestica
A04     Malus_domestica
A50     Outgroup
A51     Outgroup
A52     Outgroup
A53     Outgroup
A54     Outgroup
Ind8    xxx ## If you want to delete this sample
...     ...
IndN    Species_n
```

```
 cat SETS_quartets_Dmin.txt
```

> P1      P2      P3      P4      Dstatistic      Z-score p-value f4-ratio        BBAA    ABBA    BABA
> Malus_domestica Malus_domestica2        Malus_sieversii Malus_sylvestris        0.141346        4.85047 1.2317e-06     0.219508 366454  736560  554126



