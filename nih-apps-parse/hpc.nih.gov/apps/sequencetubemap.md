
















sequencetubemaps: discovering fusion genes from paired-end RNA sequencing data


  









|  |  |  |
| --- | --- | --- |
|  | 
[Biowulf High Performance Computing at the NIH](https://hpc.nih.gov)
 | 








[GitHub](https://github.com/NIH-HPC)
[YouTube](https://www.youtube.com/channel/UCx-kNd1kBskYr5KLT9-Erew)
[@nih_hpc](https://twitter.com/nih_hpc)
[RSS Feed](/hpc_RSS.xml)
 |




* [Systems](https://hpc.nih.gov/systems/)
* [Applications](https://hpc.nih.gov/apps/)
* [Reference Data](https://hpc.nih.gov/refdb/)
* [Storage](https://hpc.nih.gov/storage/)
* [User Guides](https://hpc.nih.gov/docs/user_guides.html)
* [Training](https://hpc.nih.gov/training/)
* [User Dashboard](https://hpcnihapps.cit.nih.gov/auth/dashboard/)
* [How To](https://hpc.nih.gov/docs/how_to.html)
* [About](https://hpc.nih.gov/about/)







**sequencetubemaps: A JavaScript module for the visualization of genomic sequence graphs.**


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
[Batch job](#sbatch) 
[Swarm of jobs](#swarm) 
 |



A JavaScript module for the visualization of genomic sequence graphs. It automatically generates a "tube map"-like visualization of sequence graphs which have been created with [vg](https://github.com/vgteam/vg).



Documentation
* [sequencetubemaps GitHub Page](https://github.com/vgteam/sequenceTubeMap)


Important Notes
* Module Name: sequencetubemaps (see [the modules page](/apps/modules.html) for more information)
 * To start seqauencetubemap please start interactive session with 2 tunnel (see [tunneling](https://hpc.nih.gov/docs/tunneling/), and copy the whole directory to local.
	+ **SEQTM** sequencetubemaps installation directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive --mem=4g --gres=lscratch:10 -T -T** 
[user@cn3144 ~]$ **module load sequencetubemaps** 

```


```

[user@cn3144 ]$ **cd /data/$USER; cp -r SEQTM .** 
[user@cn3144 ]$ **cd sequenceTubeMap**
[user@cn3144 ]$ **npm start**
> sequence-tube-maps@0.1.0 start

[backend]
[backend] > sequence-tube-maps@0.1.0 serve
[backend] > node ./src/server.mjs "42583"
[backend]
[backend] TubeMapServer listening on http://:::3000
[frontend] Attempting to bind to HOST environment variable: 127.0.0.1
[frontend] If this was unintentional, check that you haven't mistakenly set it in your shell.
[frontend] Learn more here: https://cra.link/advanced-config
[frontend]
[frontend] [HPM] Proxy created: /api  -> http://localhost:3000
[frontend] [HPM] Proxy created: /api/v0  -> http://localhost:3000
[frontend] (node:1074683) [DEP_WEBPACK_DEV_SERVER_ON_AFTER_SETUP_MIDDLEWARE] DeprecationWarning: 'onAfterSetupMiddleware' option is deprecated. Please use the 'setupMiddlewares' option.
[frontend] (Use `node --trace-deprecation ...` to show where the warning was created)
[frontend] (node:1074683) [DEP_WEBPACK_DEV_SERVER_ON_BEFORE_SETUP_MIDDLEWARE] DeprecationWarning: 'onBeforeSetupMiddleware' option is deprecated. Please use the 'setupMiddlewares' option.
[frontend] Starting the development server...
[frontend]
[frontend] Compiled successfully!
[frontend]
[frontend] You can now view sequence-tube-maps in the browser.
[frontend]
[frontend]   http://127.0.0.1:41719
[frontend] Note that the development build is not optimized.
[frontend] To create a production build, use npm run build.
[frontend]
[frontend] webpack compiled successfully


```

#wait until the server start successfully, then copy the url to web browser, eg: http://127.0.0.1:41719 
![seqtm_web_browser](/images/seqtm.png)










[HPC @ NIH](https://hpc.nih.gov)  ~
[Contact](https://hpc.nih.gov/about/contact.html)


[Disclaimer](https://hpc.nih.gov/docs/disclaimer.html) ~ 
[Privacy](https://hpc.nih.gov/docs/privacy.html) ~ 
[Accessibility](https://hpc.nih.gov/docs/accessibility.html) ~ 
[CIT](https://cit.nih.gov/) ~ 
[NIH](https://www.nih.gov/) ~ 
[DHHS](https://www.dhhs.gov/) ~ 
[USA.gov](https://www.firstgov.gov/) ~
[HHS Vulnerability Disclosure](https://www.hhs.gov/vulnerability-disclosure-policy/index.html)



  



