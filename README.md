# Završni rad
**Pristupnik:** Karlo Jularić (0036558556)

**Studij:** Elektrotehnika i informacijska tehnologija i Računarstvo

**Modul:** Računarstvo

**Mentor:** doc. dr. sc. Krešimir Križanović

**Naslov teme:** Analiza metagenomskog uzorka na temelju distribucije minimizatora

## Opis teme
U bioinformatici, metagenomski uzorak sadrži genetski materijal više različitih organizama. Osnovni korak u analizi metagenomskog uzorka je određivanje organizama koji se u njemu nalaze te njihove zastupljenosti u uzorku. U analizi nizova, k-mer je podniz duljine k. Minimizator je abecedno najmanji k-mer u kraćoj podskevenci. Za odabranu duljinu k, za svaku sekvencu je moguće odrediti distribuciju svih k-mera kao i distribuciju svih minimizatora koji se u toj sekvenci pojavljuju.
U sklopu ovog rada potrebno je ispitati mogućnost upotrebe distribucije minimizatora za analizu metagenomskog uzorka. Koristiti nekoliko genoma iz javno dostupne baze podataka RefSeq. Simulirati očitanja te generirati metagenomski uzorak. Očitanja iz metagenomskog uzorka klasificirati na temelju sličnosti distribucije minimizatora između očitanja i genoma.
Programski kod je potrebno komentirati i pri pisanju pratiti neki od standardnih stilova. Napisati iscrpne upute za instalaciju i korištenje. Kompletno programsko rješenje postaviti na GitHub.

## Zavisnosti

### C++
- CMake 3.15+
- C++17 kompatibilni prevoditelj (g++ ili clang++)

### Python
- Python 3.8+
- pandas
- scikit-learn

### R
- R 4.0+
- ggplot2
- reshape2
- caret
- rmarkdown

## Instalacija

### Python okruženje
```bash
python3 -m venv venv
source venv/bin/activate
pip install pandas scikit-learn
```

### R paketi
```r
install.packages(c("ggplot2", "reshape2", "caret", "rmarkdown"))
```

### Kompajliranje C++ programa
```bash
mkdir build && cd build
cmake ..
make
```

## Konfiguracija

Program se konfigurira putem datoteke `config.txt` u korijenskom direktoriju projekta:

```
k=6
w=10
reference_files=genome1.fasta,genome2.fasta
fragment_files=reads1_trimmed.fastq,reads2_trimmed.fastq
```

- `k` — duljina k-mera
- `w` — veličina prozora za odabir minimizatora
- `reference_files` — nazivi FASTA datoteka referentnih genoma (smještenih u `data/reference_files/`)
- `fragment_files` — nazivi FASTQ datoteka fragmenata (smještenih u `data/fragments_files/`)

## Pokretanje

```bash
cd build
./distribution_analysis
```

Program automatski pokreće Python skriptu za evaluaciju i R skriptu za generiranje PDF izvještaja.

## Output

Nakon izvršavanja program generira sljedeće datoteke u direktoriju `output/`:

| Datoteka | Opis |
|---|---|
| `reference_data.csv` | Distribucijski vektori referentnih genoma |
| `fragment_data.csv` | Distribucijski vektori fragmenata |
| `classification_data.csv` | Rezultati klasifikacije — fragment, najsličniji referentni genom, stvarna referenca i kosinusna sličnost |
| `confusion_matrix.csv` | Matrica zabune generirana Python skriptom |
| `metagenomic_reference.fasta` | Spojeni referentni genomi |
| `distribution_report.pdf` | PDF izvještaj s vizualizacijama |
