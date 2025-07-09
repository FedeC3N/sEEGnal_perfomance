### sEEGnal

Package for automatic preprocessing of EEG recordings.

This repository comprises three high-level blocks, namely, standardize, badchannel detection and artifact detection.

## Standardize

Converts the original data structure to a BIDS-compliance data structure.

For more information regarding BIDS, consult [BIDS official website](https://bids-specification.readthedocs.io/en/stable/index.html).

## Badchannel detection

Marks badchannels in EEG recordings based on different criteria.

For EEG:

- Channels with impedances above a certain threshold.
- Channels with impossible low or high amplitudes.
- Channels with significantly higher energy in 45-55 Hz range compared to the rest of channels.
- Channels with gel bridge.
- Channels with significantly higher standard deviation of amplitude compared to the rest of channels.

## Artifact detection

First, performs an independent component analysis (ICA) and then label the ICs using [MNE-ICALabel](https://mne.tools/mne-icalabel/stable/index.html).

Then looks for EOG arrtifacts, muscle artifacts, sensor artifacts, and "other" artifacts.

- EOG. Filter the recording in low frequencies. Compare frontal channels vs rest of the channels and look for high amplitude peaks significantly different.
- Muscle. Get the time series reconstructed using "muscle" components. Filter between 110-145 Hz and look for high amplitude bursts.
- Sensor. Filter the data in low frequencies. Look for high amplitude peaks.
- Other. Look for bursts with impossible amplitudes.


## Acknoledgemnents

This repository has been developed under the European project AI-Mind (European Union’s Horizon 2020 research and innovation programme, grant agreement No 964220).

## Contributors
[<img src="https://avatars.githubusercontent.com/u/138225612?size=80">](https://github.com/FedeC3N)     [<img src="https://meg.ucm.es/wp-content/uploads/2019/10/Ricardo-Bru%C3%B1a-e1592997594200-225x300.jpg? " width="80" height="80">](https://github.com/rbruna)     [<img src="https://www.gravatar.com/avatar/a49420a87cbf659b27a78c921b94cc8b?s=80&d=identicon">](https://gitlab.lurtis.com/v.ayllon)