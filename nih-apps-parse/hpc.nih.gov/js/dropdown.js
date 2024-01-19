( function( $ ) {
$( document ).ready(function() {
// Cache the elements we'll need
var menu = $('#cssmenu');
var menuList = menu.find('ul:first');
var listItems = menu.find('li').not('#responsive-tab');
//var baseURL = window.location.origin;
var baseURL = 'https://hpc.nih.gov';

// Create responsive trigger
menuList.prepend('<li id="responsive-tab"><a href="#">Menu</a></li>');

// Toggle menu visibility
menu.on('click', '#responsive-tab', function(){
	listItems.slideToggle('fast');
	listItems.addClass('collapsed');
});

// This script defines the dropdown menu links
//=============================================================================

// Links under the Status header
var Systems_menulist='' +
    '<ul>' +
    '    <li><a class="sublink" href="' + baseURL + '/systems/hardware.html"			><span>Hardware</span></a></li>' +
    '    <li><a class="sublink" href="' + baseURL + '/systems/status"				><span>Status</span></a></li>' +
    '    <li><a class="sublink" href="' + baseURL + '/systems/status/index.php#part_status"	><span>Partition Status</span></a></li>' +
    '    <li><a class="sublink" href="' + baseURL + '/systems/status/partitions.html"		><span>Usage History</span></a></li>' +
    '    <li><a class="sublink" href="' + baseURL + '/systems/status/index.php#batch_limits"	><span>Batch Limits</span></a></li>' +
    '    <li><a class="sublink" href="' + baseURL + '/systems/status/index.php#licenses"	><span>Licenses</span></a></li>' +
    '</ul>';

// Links under the Applications header
var Applications_menulist='' +
    '<ul>' +
    '    <li><a class="sublink" href="' + baseURL + '/apps/#compchem"				><span>Computational Chemistry</span></a></li>' +
    '    <li><a class="sublink" href="' + baseURL + '/apps/#deeplearning"			><span>Deep Learning</span></a></li>' +
    '    <li><a class="sublink" href="' + baseURL + '/development/"				><span>Development Tools</span></a></li>' +
    '    <li><a class="sublink" href="' + baseURL + '/apps/#hts"				><span>High-Throughput Sequencing</span></a></li>' +
    '    <li><a class="sublink" href="' + baseURL + '/apps/#image"				><span>Image Analysis</span></a></li>' +
    '    <li><a class="sublink" href="' + baseURL + '/apps/#linkage"				><span>Linkage/Phylogenetics</span></a></li>' +
    '    <li><a class="sublink" href="' + baseURL + '/apps/#math"				><span>Mathematics/Statistics</span></a></li>' +
    '    <li><a class="sublink" href="' + baseURL + '/apps/#compchem"				><span>Molecular Modeling</span></a></li>' +
    '    <li><a class="sublink" href="' + baseURL + '/apps/#massspec"				><span>Mass Spectrometry</span></a></li>' +
    '    <li><a class="sublink" href="' + baseURL + '/apps/#seqanal"				><span>Sequence Analysis</span></a></li>' +
    '    <li><a class="sublink" href="' + baseURL + '/refdb/"					><span>Scientific Databases</span></a></li>' +
    '    <li><a class="sublink" href="' + baseURL + '/apps/#structbio"				><span>Structural Biology</span></a></li>' +
    '    <li><a class="sublink" href="' + baseURL + '/apps/#sysbiol"				><span>Systems Biology</span></a></li>' +
    '    <li><a class="sublink" href="' + baseURL + '/apps/#utils"				><span>Utilities</span></a></li>' +
    '</ul>';

// Links under the Storage header
var Storage_menulist='' +
    '<ul>' +
    '    <li><a class="sublink" href="' + baseURL + '/storage/index.html"			><span>Storage Overview</span></a></li>' +
    '    <li><a class="sublink" href="' + baseURL + '/storage/sharing_data.html"		><span>Sharing Data</span></a></li>' +
    '    <li><a class="sublink" href="' + baseURL + '/storage/permissions.html"			><span>Access Permissions</span></a></li>' +
    '    <li><a class="sublink" href="' + baseURL + '/storage/backups.html"			><span>Backups</span></a></li>' +
    '    <li><a class="sublink" href="https://hpcnihapps.cit.nih.gov/auth/dashboard/storage_request.php"		><span>Storage Request Form</span></a></li>' +
    '    <li><a class="sublink" href="https://hpcnihapps.cit.nih.gov/auth/dashboard/shared_data_request.php"	><span>Shared Data Request Form</span></a></li>' +
    '    <li><a class="sublink" href="https://hpcnihapps.cit.nih.gov/auth/dashboard/group_change.php"   ><span>Group Change Request Form</span></a></li>' +
//    '    <li><a class="sublink" href="' + baseURL + '/storage/object.html"			><span>Object Storage</span></a></li>' +
//    '    <li><a class="sublink" href="' + baseURL + '/nih/object_request.html"			><span>Object Storage Request Form</span></a></li>' +
    '</ul>';

// Links under the User Guides header
var User_Guides_menulist='' +
    '<ul>' +
    '    <li><a class="sublink" href="' + baseURL + '/docs/userguide.html"			><span>Biowulf User Guide</span></a></li>' +
    '    <li><a class="sublink" href="' + baseURL + '/docs/biowulf_tools.html" 			><span>Biowulf Tools</span></a></li>' +
    '    <li><a class="sublink" href="' + baseURL + '/docs/deep_learning.html"			><span>Deep Learning</span></a></li>' +
    '    <li><a class="sublink" href="' + baseURL + '/development/"				><span>Development Tools</span></a></li>' +
    '    <li><a class="sublink" href="' + baseURL + '/docs/globus/setup.php"			><span>Globus</span></a></li>' +
    '    <li><a class="sublink" href="' + baseURL + '/docs/diy_installation"			><span>Personal Software Installation</span></a></li>' +
    '    <li><a class="sublink" href="' + baseURL + '/apps/singularity.html"			><span>Singularity Containers</span></a></li>' +
    '    <li><a class="sublink" href="' + baseURL + '/docs/tunneling"				><span>SSH Tunneling</span></a></li>' +
    '    <li><a class="sublink" href="' + baseURL + '/apps/swarm.html"				><span>Swarm</span></a></li>' +
    '    <li><a class="sublink" href="' + baseURL + '/apps/jupyter.html"			><span>Jupyter Notebooks</span></a></li>' +
    '    <li><a class="sublink" href="' + baseURL + '/apps/modules.html"			><span>Environment Modules</span></a></li>' +
    '    <li><a class="sublink" href="' + baseURL + '/docs/FAQ.html"				><span>FAQ</span></a></li>' +
    '</ul>';

// Links under the Training header
var Training_menulist='' +
    '<ul>' +
    '    <li><a class="sublink" href="' + baseURL + '/training/intro_biowulf"			><span>Biowulf Online Class</span></a></li>' +
    '    <li><a class="sublink" href="' + baseURL + '/training/#seminar"			><span>Seminar Series</span></a></li>' +
    '    <li><a class="sublink" href="' + baseURL + '/training/#upcoming"			><span>Upcoming Classes</span></a></li>' +
    '    <li><a class="sublink" href="' + baseURL + '/training/#past_classes"			><span>Past Classes</span></a></li>' +
    '    <li><a class="sublink" href="' + baseURL + '/docs/trainingvids.html"			><span>Training Videos</span></a></li>' +
    '    <li><a class="sublink" href="' + baseURL + '/nih/student.html"				><span>Student Accounts (NIH only)</span></a></li>' +
    '</ul>';

// Links under the How To header
var How_To_menulist='' +
    '<ul>' +
    '    <li><a class="sublink" href="' + baseURL + '/docs/accounts.html"			><span>Get An Account</span></a></li>' +
    '    <li><a class="sublink" href="' + baseURL + '/docs/how_to.html"				><span>Unlock Your Acccount</span></a></li>' +
    '    <li><a class="sublink" href="' + baseURL + '/docs/connect.html"			><span>Connect To Helix/Biowulf</span></a></li>' +
    '    <li><a class="sublink" href="' + baseURL + '/docs/transfer.html"			><span>Transfer Files</span></a></li>' +
    '    <li><a class="sublink" href="https://hpcnihapps.cit.nih.gov/auth/dashboard/speedtest/"			><span>Estimate Transfer Speeds</span></a></li>' +
    '    <li><a class="sublink" href="https://hpcnihapps.cit.nih.gov/auth/dashboard/storage_request.php"		><span>Request More Disk Space</span></a></li>' +
    '    <li><a class="sublink" href="' + baseURL + '/storage/sharing_data.html"		><span>Share Data</span></a></li>' +
    '    <li><a class="sublink" href="https://hpcnihapps.cit.nih.gov/auth/dashboard/group_change.php"   ><span>Group Change Request Form</span></a></li>' +
    '    <li><a class="sublink" href="' + baseURL + '/docs/how_to.html"				><span>and more...</span></a></li>' +
    '</ul>';

// Links under the About header
var About_menulist='' +
    '<ul>' +
    '    <li><a class="sublink" href="' + baseURL + '/systems/"					><span>Systems</span></a></li>' +
    '    <li><a class="sublink" href="' + baseURL + '/about/biowulf_stats.html"			><span>Biowulf Stats</span></a></li>' +
    '    <li><a class="sublink" href="' + baseURL + '/about/index.html#staff"			><span>HPC Staff</span></a></li>' +
    '    <li><a class="sublink" href="' + baseURL + '/docs/userguide.html#ack"			><span>Citation</span></a></li>' +
    '    <li><a class="sublink" href="' + baseURL + '/about/history.html"			><span>History</span></a></li>' +
    '    <li><a class="sublink" href="' + baseURL + '/about/recruit.html"			><span>Careers</span></a></li>' +
    '    <li><a class="sublink" href="' + baseURL + '/policies/index.html"			><span>Policies</span></a></li>' +
    '    <li><a class="sublink" href="' + baseURL + '/about/contact.html"			><span>Contact</span></a></li>' +
    '</ul>';

// Now append the contents above into the page
$('#Systems_menulist').append(Systems_menulist);
$('#Applications_menulist').append(Applications_menulist);
$('#Storage_menulist').append(Storage_menulist);
$('#User_Guides_menulist').append(User_Guides_menulist);
$('#Training_menulist').append(Training_menulist);
$('#How_To_menulist').append(How_To_menulist);
$('#About_menulist').append(About_menulist);

//=============================================================================

});




} )( jQuery );

