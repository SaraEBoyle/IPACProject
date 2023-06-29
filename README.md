# IPACProject
Contains MATLAB, Arduino, and Bonsai code used to collect and analyze data published in "Neurotensin neurons in the extended amygdala control dietary choice and energy homeostasis," 2022 from Furlan et. al.

Analysis package for the 2-spout delivery task, a task designed together with Alessandro Furlan and implemented by Sara Boyle. In this task, 2 different liquids are loaded and calibrated for delivery to a mouse. When a trial starts, a droplet of one of the liquids is randomly dispensed to sit at the end of a spout. The mouse can then freely approach and drink the fluid. When the mouse licks the spouts, it closes a simple circuit, allowing the Bpod to register the timing of the individual spout licks, and triggering the next trial. This analysis package uses the timing data from the Bpod and photometry data from the Neurofiberphotometrics photometry system to analyze photometry responses from mice to different types and volumes of liquid stimuli. It includes built in bleaching correction and many customizable parameters.

<img width="65" alt="Screen Shot 2023-06-28 at 12 35 03 PM" src="https://github.com/SaraEBoyle/IPACProject/assets/83416542/3c05ae9c-8985-4250-9db6-e19340843a0a">

Inputs (must all be contained in same folder):
--> photometry data (.csv), collected through custom Bonsai program, using Neurofiberphotometrics photometry system.
--> analog timing information (.csv), collected through custom Bonsai program. Contains analog input recorded from Bpod system, marking the start of liquid delivery trials. Needed to synchronize photometry data to behavioral data.
--> experimental behavioral data (.mat) from Bpod 2 spout protocol. Includes trial start timing information and timing on when mice lick the delivery spouts.

Outputs:
<-- Segments the channels recorded by the photometry system and corrects the session for bleaching effects
<-- For each type of liquid and delivery volume specified in Bpod protocol, segments photometry data and plots heatmaps, average dF/F, and Z-Scores of dF/F.
<-- Also calculates licking rate for each liquid delivered.

1. Run Two_Spout_Analysis.m.
2. The program gives the option to record photometry data (with patch cords attached) without the system LEDs turned on in order to subtract ambient light. If you choose to do this (recommended), follow the prompts to select the period of data to use for subtractions.

![raw photometry data](https://github.com/SaraEBoyle/IPACProject/assets/83416542/f0c09387-f111-4b9e-bed0-61ef6db973cc)

3. Afterwards, you will have a chance to trim the data to cut away anything from the beginning or end of your session that may influence bleaching correction.

![raw data after cutting beginning](https://github.com/SaraEBoyle/IPACProject/assets/83416542/cac872e2-965b-4c9c-ac51-64b5a835236d)

5. Once the data is appropriately trimmed, it will be split into the correct channels. This code is built to handle 2 channels, a control, isosbestic channel, and a signal channel.

   ![Calculate data split](https://github.com/SaraEBoyle/IPACProject/assets/83416542/f15a4c44-a7fb-4be1-a72c-f4a9a9ddd448)

![split channels](https://github.com/SaraEBoyle/IPACProject/assets/83416542/6c799ff5-b5b4-418c-b5e0-776ab495ef2e)

5. Then, it will apply bleaching correction to both channels. You have the option to use a biexponential curve fit (recommended) or a sliding window method.

![bleaching correction fit](https://github.com/SaraEBoyle/IPACProject/assets/83416542/8fb1ef6b-7fb1-4719-9099-414fdd32401a)
![after bleaching correction](https://github.com/SaraEBoyle/IPACProject/assets/83416542/21eea6bf-ba35-4045-8dc9-735199a0ca26)

6. Once both channels are corrected, it will save the signal and control channels and display the entire session for you.

![both channels after bleaching correction](https://github.com/SaraEBoyle/IPACProject/assets/83416542/7b36b758-57c1-451b-957b-4938ed9e1638)

7. From there, the program will use the timing data from the Bpod file to align the photometry data to relevant liquid delivery events, plotting heatmaps, average traces, and lick rate informrmation from all trial types.

![photometry response heatmap](https://github.com/SaraEBoyle/IPACProject/assets/83416542/94f80e82-b3a5-41b2-8c4f-a6aceae3d11f)

![photometry response average](https://github.com/SaraEBoyle/IPACProject/assets/83416542/d3c8ebe9-ca61-45fc-bf96-d7446ca5750f)



