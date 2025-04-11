### AI-Mind ETL process 

![Python Logo](https://www.ai-mind.eu/wp-content/uploads/sites/39/2020/09/Al-Mind_logo_web-300x185.png "Sample inline image")

Package for automatic preprocessing of EEG and MEG recordings.

This repository comprises four high-level blocks, namely, standardize, badchannel detection, artifact detection, and quality assesment.

## Standardize

Converts the original data structure to a BIDS-compliance data structure.

For more information regarding BIDS, consult [BIDS official website](https://bids-specification.readthedocs.io/en/stable/index.html).

## Badchannel detection

Marks badchannels in EEG and MEG recordings based on different criteria.

For EEG:

- Channels with impedances above 200 KΩ.
- Channels with significantly higher variance compared to the rest of channels.
- Channels with significantly higher energy in 45-55 Hz range compared to the rest of channels.
- Channels with gel bridge.

For MEG:

- Channels with no amplitude due to tunnings problems (flat channels).
- Channels with significantly higher variance compared to the rest of channels.
- Channels with significantly higher energy in 45-55 Hz range compared to the rest of channels.

## Artifact detection

First, performs an independent component analysis (ICA) and then label the ICs using [MNE-ICALabel](https://mne.tools/mne-icalabel/stable/index.html).

Removes the components labelled as EOG, EKG, muscle, line noise, and channel noise.

Looks for remaining artifacts, specifically, it looks for EOG arrtifacts, muscle artifacts, and sensor artifacts (jumps).

## Quality assessment

Checks the EEG/MEG recordings quality after processing based on:

- if any files has been corrupted through the process.
- the number of badchannels.
- the number of clean epochs.

If the recording fulfill the quailty assessment, it is stored in a different folder with BIDS structure.

## Acknoledgemnents

This repository has been developed under the European projectg AI-Mind (European Union’s Horizon 2020 research and innovation programme, grant agreement No 964220).

## Contributors
[<img src="https://avatars.githubusercontent.com/u/138225612?size=80">](https://github.com/FedeC3N)     [<img src="https://meg.ucm.es/wp-content/uploads/2019/10/Ricardo-Bru%C3%B1a-e1592997594200-225x300.jpg? " width="80" height="80">](https://github.com/rbruna)     [<img src="https://www.gravatar.com/avatar/a49420a87cbf659b27a78c921b94cc8b?s=80&d=identicon">](https://gitlab.lurtis.com/v.ayllon)