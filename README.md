### AI-Mind ETL process 

![Python Logo](https://www.ai-mind.eu/wp-content/uploads/sites/39/2020/09/Al-Mind_logo_web-300x185.png "Sample inline image")

Package for automatic preprocessing of EEG recordings.

This repository comprises four high-level blocks, namely, standardize, badchannel detection, and artifact detection.

## Standardize

Converts the original data structure to a BIDS-compliance data structure.

For more information regarding BIDS, consult [BIDS official website](https://bids-specification.readthedocs.io/en/stable/index.html).

## Badchannel detection

Marks badchannels in EEG recordings based on different criteria:

- Channels with impedances above 200 KΩ.
- Channels with significantly higher variance in amplitude compared to the rest of channels.
- Channels with significantly higher energy in 45-55 Hz range compared to the rest of channels.
- Channels with gel bridge.


## Artifact detection

Since muscle and sensor artifacts may lead the estimation of Ics, we perform two consecutive SOBIs:
-	The first SOBI (and the associated ICLabel) is used to identify epochs with muscular and sensor artifacts. Before estimation SOBI for the second time, these muscular and sensor artifacts are removed.
-	The second SOBI (and the associated ICLabel) is used to identify the remaining muscle, sensor, and ocular artifacts.



## Acknoledgemnents

This repository has been developed under the European projectg AI-Mind (European Union’s Horizon 2020 research and innovation programme, grant agreement No 964220).

## Contributors
[<img src="https://avatars.githubusercontent.com/u/138225612?size=80">](https://github.com/FedeC3N)     [<img src="https://meg.ucm.es/wp-content/uploads/2019/10/Ricardo-Bru%C3%B1a-e1592997594200-225x300.jpg? " width="80" height="80">](https://github.com/rbruna)     [<img src="https://www.gravatar.com/avatar/a49420a87cbf659b27a78c921b94cc8b?s=80&d=identicon">](https://gitlab.lurtis.com/v.ayllon)