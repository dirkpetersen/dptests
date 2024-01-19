

document.querySelector('title').textContent = ' aria2 on Biowulf ';

aria2 on Biowulf 


Description

`aria2c` is a utility for downloading files over the HTTP(S), (S)FTP,
and Metalink protocols. The version installed on our systems does not support
bittorrent. In can use multiple sources and multiple connections to the same
source to speed up downloads. It can make use of the proxies that allow limited
internet access to the compute nodes.



There may be multiple versions of aria2 available. An easy way of selecting the
version is to use [modules](/apps/modules.html). To see the modules
available, type



```

module avail aria2 

```

To select a module use



```

module load aria2/[version]

```

where `[version]` is the version of choice.


### Environment variables set


* $PATH


### Documentation


* [Manual](https://aria2.github.io/manual/en/html/index.html)
* [Github](https://github.com/aria2/aria2)




Usage
Download a file with a single connection from a single source



```

$ module load aria2
$ aria2c "$testurl"

```

Download a file with 4 single connections from a single source. `-x`
is used to set the number of connections and `-k` is used to
specify the chunk size.



```

$ aria2c -x4 -k1M "$testurl"

```

Download files listed in a text file (`-i`) with a max of 2 concurrent 
downloads (`-j`).



```

$ aria2c -i listfile -j2

```

Throttle download speed per download



```

$ aria2c --max-download-limit=1M $testurl

```

Throttle overall download speed



```

$ aria2c --max-overall-download-limit=1M $testurl

```





