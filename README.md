# AnalyzeIUPAC
short python program to analyze some iupac names     
\* not intended for serious use, just made this really quick for fun

# Usage
```python
>> from analyzeiupac import Molecule
>> name = "6,6-dimethyl-4-methylenebicyclo[3.1.1]hept-2-en-1-ol"
>> data = Molecule(name).data
parent: hept
cyclo: bicyclo[3.1.1]
unsaturation: ['2', 'en']
functionalgroup: ['1', 'ol']
substituents: 6,6-dimethyl-4-methylene
>> print(data)
('hept', 'bicyclo[3.1.1]', ['2', 'en'], ['1', 'ol'], '6,6-dimethyl-4-methylene')
```

# Info
- put the location of the unsaturation directly infront of the "en"/"yn"/"dien"/"diyn"/etc. (as shown in the example above)     
- cyclo, unsaturation, functional group, and substituents are all optional. The program uses the position of the parent string (and cyclo, if found) to analyze the rest of the molecule.
- trans/cis is ignored
