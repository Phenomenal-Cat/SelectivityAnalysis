# Feature Selectivity Analysis
<b>Matlab code for analysis and graphical representation of single-unit neurophysiological data collected during visual presentation of parameterized 3D stimuli.</b>

For social animals like primates, faces are a behaviourally important class of visual objects, and they exhibit statistical variation across multiple dimensions. Previous studies have typically investigated the selectivity of single neurons in the inferotemporal (IT) cortex of macaques to variation of one or two facial features at a time, such as identity and head orientation. However, due to the high-dimensionality of the stimulus space and the time-constraints imposed by traditional recording methods, our understanding of how multiple facial features are encoded remains limited. 

<img src="http://jn.physiology.org/content/jn/112/7/1748/F1.large.jpg" width=400 height=300> 
<b>Figure 1.</b> (A-D) Chronic microwire elctrode, figure courtesy of McMahon et al. (2014). 

To overcome these limitations and obtain a broader characterization of IT neuronal preferences for visual face features, we exploit two recent methodological advancements. First, we surgically implanted chronic microwire ‘brush’ electrodes, which allow us to stably isolate single units over periods of weeks to months. The ability to record from the same neuron longitudinally increases the total number of visual stimuli that can be tested, thus improving coverage of stimulus space and characterization of the cell’s selectivity across multiple stimulus dimensions. To leverage this, we have developed an anatomically realistic three-dimensional virtual macaque avatar that offers continuous parametric control of many variables (in addition to lighting, material and environmental varibles). To sample this high-dimensional feature space we use an adaptive approach whereby neuronal responses from one session are used to inform which regions of feature space to explore in the next session, with each generation of new stimuli being rendered offline using ray-tracing. 

<img src="https://cloud.githubusercontent.com/assets/7523776/25898502/5c17dc8c-355a-11e7-91ed-9ceb096962eb.gif" width=600 height=600>
<b>Figure 2.</b> Head orientation tuning in an example AM neuron.


<b>Murphy AP & Leopold DA (2017).</b> Measuring neuronal selectivity for facial features in macaque inferotemporal cortex through adaptive sampling of feature space. Society for Neuroscience Abstract
