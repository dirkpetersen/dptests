

document.querySelector('title').textContent = 'Text Editors';
php
 /\* function for obtaining version numbers, David 10/29/09 \*/
 include('print\_version.php');
?
Text Editors
A variety of text and source code editors are available for both X-windows and text terminals.


**Nano (Pico)**


Pico and its clone, nano, are simple, user-friendly text
editors derived from the editor in the Pine email client. Type *pico [filename]* or *nano [filename]* to
edit a file and use the key commands listed at the bottom of the screen to
access various functions.   


[[Summary of Pico commands](http://www.csee.umbc.edu/courses/104/spring05/dblock/picocmds.shtml)]
[[Pico command-line options](https://en.wikipedia.org/wiki/Pico_%28text_editor%29#Command-line_options)]




**Atom** 


Atom is a modern desktop text editor. The stated
 goal is a zero-compromise combination of hackability and usability: an
 editor that will be welcoming to an elementary school student on their
 first day learning to code, but also a tool they won't outgrow as they
 develop into seasoned hackers. Atom is a specialized variant of
 Chromium designed to be a text editor rather than a web browser. Every
 Atom window is essentially a locally-rendered web page.


To run Atom, you must use a graphical connection to Biowulf, set up an [interactive session](https://hpc.nih.gov/docs/userguide.html#int), and load the atom [module](https://hpc.nih.gov/apps/modules.html) before executing the atom command.

When you launch Atom for the first time, you should
 get a screen presenting several options for documentation and for
 customization and operation of the editor. Atom's primary
 documentation is the [[Atom
 Flight Manual](http://flight-manual.atom.io/)], including its [[Atom
 Basics](http://flight-manual.atom.io/getting-started/sections/atom-basics)] material.


**Emacs** 


Emacs is a text and source code editor for text
terminals and X. It has a vast set of features and is well suited for doing
everything from reading mail and simple text editing to managing and
editing large programming projects. It has its own help and tutorial which
can be accessed by typing *Ctrl-h i* and *Ctrl-h t* respectively. Type *emacs
[filename]* to edit a file. [[Emacs
Manual](http://www.gnu.org/software/emacs/manual/html_mono/emacs.html)]


Newer versions of emacs are available through the environment modules. See
*module avail emacs* for mor details.


**ESS**


Emacs Speaks Statistics is an Emacs mode for interactive statistical programming and data analysis. Languages supported: the S family (S, S-PLUS and R), SAS, BUGS/JAGS, Stata and XLispStat. First load an R module per [our R applications page](http://helix.nih.gov/Applications/R.html).
 Putting the line *(load "/usr/local/share/emacs/site-lisp/ess-17.11/lisp/ess-site")* in your .emacs file, or its one-time equivalent *M-x load-library /usr/local/share/emacs/site-lisp/ess-17.11/lisp/ess-site*, will make an \*ESS\* buffer available. [[ESS Documentation](http://ess.r-project.org/index.php?Section=documentation&subSection=manuals)]


**Nedit** 


NEdit is an GUI style editor for plain text and source code files. It provides mouse based editing and a streamlined editing style, based on popular Macintosh and MS Windows editors, using the X-window system. NEdit requires an X-based workstation or X-Terminal. Type *nedit [filename]* to edit a file. 
[[Nedit manual](ftp://ftp.fu-berlin.de/unix/editors/nedit/contrib/misc/nedit.pdf)]



**SciTE**


[SciTE](http://www.scintilla.org/SciTE.html) is a source code editor for X-windows that is well suited for scripting and programming small projects. Type *SciTE [filename]* to edit a file. 
[[SciTE documentation](http://www.scintilla.org/SciTEDoc.html)]



**vi**


The venerable vi editor is a text and source code editor for text terminals and X-windows. It has a large set of features for simple text editing and for sophisticated programming projects. Type *vi [filename]* to edit a file. Type *:help* to start learning vi. 
[[vi documentation](http://tuxgraphics.org/~guido/vi/)]



**vim**


vim is a text editor that is upwards compatible to Vi. It can be used to edit all kinds of plain text. It is especially useful for editing programs with syntactical coloring. There are a lot of enhancements above Vi: multi level undo, multi win- dows and buffers, syntax highlighting, command line editing, filename completion, on-line help, visual selection, etc.
[[vim documentation](http://www.vim.org/docs.php)]


**[VScode](vscode.html)**


Visual Studio Code is a lightweight but powerful source code editor which runs on your desktop and is available for Windows, macOS and Linux. It comes with built-in support for JavaScript, TypeScript and Node.js and has a rich ecosystem of extensions for other languages (such as C++, C#, Java, Python, PHP, Go) and runtimes (such as .NET and Unity). 
[[VS Code documentation](https://code.visualstudio.com/docs)]






