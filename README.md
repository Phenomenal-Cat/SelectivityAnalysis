# Feature Selectivity Analysis
<b>Matlab code for analysis and graphical representation of single-unit neurophysiological data collected during visual presentation of parameterized 3D stimuli.</b>

For social animals like primates, faces are a behaviourally important class of visual objects, and they exhibit statistical variation across multiple dimensions. Previous studies have typically investigated the selectivity of single neurons in the inferotemporal (IT) cortex of macaques to variation of one or two facial features at a time, such as identity and head orientation. However, due to the high-dimensionality of the stimulus space and the time-constraints imposed by traditional recording methods, our understanding of how multiple facial features are encoded remains limited. 

<dl><img src="https://user-images.githubusercontent.com/7523776/41604948-608d24a0-73ae-11e8-80d9-72db5f19a0d5.png" height=260> 
</br><b>Figure 1.</b> Chronic microwire elctrode location and 'brush' tip illustration (MR image courtesy of McMahon et al., 2014). 
</dl>

To overcome these limitations and obtain a broader characterization of IT neuronal preferences for visual face features, we exploit two recent methodological advancements. First, we surgically implanted chronic microwire ‘brush’ electrodes, which allow us to stably isolate single units over periods of weeks to months. The ability to record from the same neuron longitudinally increases the total number of visual stimuli that can be tested, thus improving coverage of stimulus space and characterization of the cell’s selectivity across multiple stimulus dimensions. 

<dl><img src="https://user-images.githubusercontent.com/7523776/41604954-643a99fc-73ae-11e8-9479-cffe8482a8d3.png" height=300>
</br><b>Figure 2.</b> "Finger printing" an example face-selective neuron across days. Each day a reduced stimulus set consisting of 60 images (10 exemplars from 6 categories) is presented at the start of the session. Even neighbouring neurons show highly stereotypical tuning across the stimuli, and the stability of this image-tuning signature can be used to ensure that we are indeed isolating the same cells across sessions. 
</dl>

Second, to leverage this new ability we have developed an anatomically realistic [three-dimensional virtual macaque avatar](https://github.com/MonkeyGone2Heaven/MacaqueBlender/wiki) that offers continuous parametric control of many variables (in addition to lighting, material and environmental varibles). To sample this high-dimensional feature space we use an adaptive approach whereby neuronal responses from one session are used to inform which regions of feature space to explore in the next session, with each generation of new stimuli being rendered offline using ray-tracing. 

<dl><img src="https://cloud.githubusercontent.com/assets/7523776/25898502/5c17dc8c-355a-11e7-91ed-9ceb096962eb.gif" width=600 height=600>
</br><b>Figure 3.</b> Head orientation tuning in an example AM neuron. This animation shows how the spike rate of a single neuron evolves over time (show in milliseconds) from when the visual stimulus appeared on the screen in front of the subject. This cell shows tuning to head orientation.
</dl>

Preliminary data collected using these methods have been reported in a scientific poster presentation, here: 

<b>Murphy AP & Leopold DA (2019).</b> [A parameterized digital 3D model of the Rhesus macaque face for investigating the visual processing of social cues. J.Neurosci.Methods. 324, 108309](https://doi.org/10.1016/j.jneumeth.2019.06.001). 

<b>Murphy AP & Leopold DA (2017).</b> [Measuring neuronal selectivity for facial features in macaque inferotemporal cortex through adaptive sampling of feature space. Society for Neuroscience Abstract](https://www.researchgate.net/publication/323126846_Measuring_neuronal_selectivity_for_facial_features_in_macaque_inferotemporal_cortex_through_adaptive_sampling_of_feature_space)
