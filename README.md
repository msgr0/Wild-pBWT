`pbwt - wildcards`
---
# Pbwt

# Blocks

# Haplotype Gen
Generate random 0-1 data in two output files: one with gaps (*) and one ground-truth file (with no gaps).

# TODO
- bitvector per ogni carattere dell'alfabeto
- durante la computazione dei blocchi, prendo il carattere successivo rispetto alla colonna dei prefissi ordinati
    - se ha un valore, cerco nel rispettivo bitvector
    - altrimenti, DEVO cercare in tutti i bit vector il fatto che sia "coerente/aperto" per almeno uno dei SIGMA-bitvector
    - se per tutti i bit vector ho "differenza" \not \equal 0 allora il blocco chiude.
    


- aggiungere PARAMETRO grandezza blocco

----
AUDIO
