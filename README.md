# pitch-detection

this code is with 3 methods
  <li> AMDF </li>
  <li> Autocorrelation </li>
  <li> Cepstrum Analysis </li>
<br>
<br>
In the self-correlation function method, 2 peaks (numerically large) of the self-correlation function diagram are selected as lag for the estimated peak.  
Figure 1: Smoothed curve of points obtained) by Center Clipping method  
Figure 2: Candidate points obtained by Center Clipping method  
Figure 3: smoothed curve of points obtained (Level Center Clipping method 3)  
Figure 4: Candidate points obtained by 3 Level Center Clipping method  
Coefficient c for Center Clipping: 0.2  
Coefficient c for 3 Level Center Clipping: 0.2  
Energy threshold for sound frame detection: 500,000  

<h2> first Method: AMDF </h2>

![image](https://user-images.githubusercontent.com/54143711/127675568-34f9aabd-408c-45e9-8704-ce98d27ebf5c.png)


<h2> second Method: Autocorrelation </h2>

![image](https://user-images.githubusercontent.com/54143711/127675855-17baa104-aa02-4fdf-b3c9-f443f152ae43.png)



<h2> first Method: Cepstrum Analysis </h2>

![image](https://user-images.githubusercontent.com/54143711/127676016-16f0170b-41d5-4db7-8d52-e298a52ac377.png)


