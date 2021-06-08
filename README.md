# Enhanced Cuckoo Search - Node Localisation (ECS-NL)

This codebase is the counterpart of the ECS-NL paper. You can read it [here](https://www.mdpi.com/1424-8220/21/11/3576).

Node localization is an essential process for setting up a Wireless Sensor Network (WSN). <br/>
This project is an application of a bio-inspired algorithm called Cuckoo Search (CS) for localizing unknown nodes in the network. <br/>

Abstract:<br/>
>Node localisation plays a critical role in setting up Wireless Sensor Networks (WSNs). A sensor in WSNs senses, processes and transmits the sensed information simultaneously. Along with the sensed information, it is crucial to have the positional information associated with the information source. A promising method to localise these randomly deployed sensors is to use bio-inspired meta-heuristic algorithms. In this way, a node localisation problem is converted to an optimisation problem. Afterwards, the optimisation problem is solved for an optimal solution by minimising the errors. Various bio-inspired algorithms, including the conventional Cuckoo Search (CS) and modified CS algorithm, have already been explored. However, these algorithms demand a predetermined number of iterations to reach the optimal solution, even when not required. In this way, they unnecessarily exploit the limited resources of the sensors resulting in a slow search process. This paper proposes an Enhanced Cuckoo Search (ECS) algorithm to minimise the Average Localisation Error (ALE) and the time taken to localise an unknown node. In this algorithm, we have implemented an Early Stopping (ES) mechanism, which improves the search process significantly by exiting the search loop whenever the optimal solution is reached. Further, we have evaluated the ECS algorithm and compared it with the modified CS algorithm. While doing so, note that the proposed algorithm localised all the localisable nodes in the network with an ALE of 0.5–0.8 m. In addition, the proposed algorithm also shows an 80% decrease in the average time taken to localise all the localisable nodes. Consequently, the performance of the proposed ECS algorithm makes it desirable to implement in practical scenarios for node localisation.

Note:
These are various versions of the ECS algorithm (Written in Python). <br/>
<b>You're looking for "ecs-basecode.ipynb".</b>
