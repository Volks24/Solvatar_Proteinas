# Solvatar_Proteinas

# generar_proteina_solvatada.py

Este script permite insertar una proteína en el centro de una caja de solventes (agua, etanol, fenol, etc.) y eliminar automáticamente aquellos solventes que se solapen con la proteína, generando así un sistema listo para simulaciones de dinámica molecular.

---

## ¿Qué hace este script?

1. **Carga una proteína en formato PDB**.
2. **Carga una caja de solventes (puede ser cúbica o octaédrica)**.
3. **Centra la proteína en el centro geométrico de la caja**.
4. **Elimina solventes que se encuentren a menos de cierta distancia de cualquier átomo de la proteína**.
5. **Genera un archivo `.pdb` con la proteína correctamente centrada y rodeada por solventes no solapados**.

---

## Requisitos

- Python 3
- [MDTraj](http://mdtraj.org/)
- NumPy
- SciPy

Podés instalar las dependencias con:

```bash
pip install mdtraj numpy scipy
```

---

## Uso

```bash
python generar_proteina_solvatada.py -p proteina.pdb -b caja_solv.pdb -d 4 -s WAT ETH PHE
```

### Parámetros:

- `-p`, `--proteina`: archivo `.pdb` de la proteína (sin agua ni ligandos).
- `-b`, `--box`: archivo `.pdb` de la caja de solventes.
- `-d`, `--distancia`: distancia mínima en Å para eliminar solventes cercanos a la proteína. *(default = 4 Å)*
- `-s`, `--solventes`: lista de `resname` de solventes a considerar para eliminar. *(default = WAT ETH PHE)*

---

## Salidas

- Un archivo llamado `<nombre_proteina>_center.pdb`: contiene la proteína centrada en el origen.
- Un archivo llamado `<nombre_proteina>_SV.pdb`: contiene el sistema final con la proteína centrada en la caja y los solventes filtrados.

---

## Autor

Juan Manuel Prieto – [@jmprietobio](https://github.com/Volks24)
