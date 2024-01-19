

document.querySelector('title').textContent = 'Stow on Biowulf';
Stow on Biowulf


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Example Session](#int) 
 |



GNU Stow is a symlink farm manager which takes distinct packages of software and/or data located in separate directories on the filesystem, and makes them appear to be installed in the same place. For example, /usr/local/bin could contain symlinks to files within /usr/local/stow/emacs/bin, /usr/local/stow/perl/bin etc., and likewise recursively for any other subdirectories such as .../share, .../man, and so on.

This is particularly useful for keeping track of system-wide and per-user installations of software built from source, but can also facilitate a more controlled approach to management of configuration files in the user's home directory, especially when coupled with version control systems.

Stow is safe-- it does not overwrite any files it didn't create.



Documentation
* [GNU Stow Main Site](https://www.gnu.org/software/stow/)
* [GNU Stow Documentation](https://www.gnu.org/software/stow/manual/)


Important Notes
* Module Name: stow (see [the modules page](/apps/modules.html) for more information)



Example Session

One good use case for stow is where you've built some software and want to symlink it to a directory that's already on your PATH.

  
Sample session (user input in **bold**):



```

[user@helix]$ **ls myapp**
bin
src
lib
include
README
[user@helix ~]$ **which myapp # not curently on the PATH**
/usr/bin/which: no myapp in (...)
[user@cn3144 ~]$ **module load stow**
[+] Loading stow, version 2.2.2...
[user@helix ~]$ **stow --target $HOME/.local myapp**
[user@helix ~]$ **which myapp**
/home/user/.local/bin/myapp
[user@helix ~]$ **stow -D --target $HOME/.local myapp # now unstow the application**
[user@helix ~]$ **which myapp # now we're back to our original state**
/usr/bin/which: no myapp in (...)
[user@helix ~]$ **logout**

```









