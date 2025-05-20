# Computational Analysis of Greek Folk Music: The case of Aegean Syrtos and Balos

This work contains a dataset of Greek Aegean folk music tunes focusing on two prominent dances, Syrtos and Balos, and its feature-pattern analysis.

## About the Dataset
The dataset consists of musical recordings and corresponding annotations for Syrtos and Balos tunes. These annotations include information about melodic and rhythmic patterns, as well as regional variations.

## Methodology
To analyze the dataset, we employed a tailored version of [PatMinr](https://github.com/olivierlar/miningsuite), a computational tool for pattern detection in music. This allowed us to identify and study musical patterns within the tunes.

## Results
Our analysis revealed distinct melodic and rhythmic variances between the Syrtos and Balos dance tune categories. Additionally, we observed clear regional influences on the musical characteristics of these dances.

## Usage
Researchers and enthusiasts interested in Greek folk music can utilize this dataset to study and explore the intricate musical nuances prevalent within folk traditions.

## Feature Analysis:

Run: MIDIfeatures.m

## Pattern Discovery Algorithm: 

in MATLAB:

Run analyzeMIDI('Syrtos-Ikaria.mid',[0,NaN,NaN,NaN,6],4) --- for patterns across one tune

Run batchanalysis --- for patterns across the corpus

You need to use the keysignature.csv!

## License
This dataset is released under the MIT license. See the LICENSE file for more details.


