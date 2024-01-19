//var host = window.location.origin;
var host = 'https://hpc.nih.gov';

// Create modification date
var months = new Array(13);
months[1] = "January";
months[2] = "February";
months[3] = "March";
months[4] = "April";
months[5] = "May";
months[6] = "June";
months[7] = "July";
months[8] = "August";
months[9] = "September";
months[10] = "October";
months[11] = "November";
months[12] = "December";
var dateObj = new Date(document.lastModified)
var lmonth = months[dateObj.getMonth() + 1]
var date = dateObj.getDate()
var fyear = dateObj.getYear()
var lhour = dateObj.getHours();
var lmin = dateObj.getMinutes();
if (fyear < 2000) 
fyear = fyear + 1900
var moddate = "Last modified: " + date + " " + lmonth + " " + fyear + " ";

// The footer
var footer = "" +
"</div>" +
"<div align='right' class='lastmod'>" + moddate + "</div>" +
"<div class='footarea'>" +
"  <div class='footer'>" +
"    <a href='" + host + "' class='footlink'>HPC @ NIH </a> ~" +
"    <a href='" + host + "/about/contact.html' class='footlink'>Contact</a>" +
"  </div>" +
"  <div class='footer2'>" +
"    <a href='" + host + "/docs/disclaimer.html' class='footlink'>Disclaimer</a> ~" +
"    <a href='" + host + "/docs/privacy.html' class='footlink'>Privacy</a> ~" +
"    <a href='" + host + "/docs/accessibility.html' class='footlink'>Accessibility</a> ~" +
"    <a href='https://cit.nih.gov/' class='footlink'>CIT</a> ~" +
"    <a href='https://www.nih.gov/' class='footlink'>NIH</a> ~" +
"    <a href='https://www.dhhs.gov/' class='footlink'>DHHS</a> ~" +
"    <a href='https://www.firstgov.gov/' class='footlink'>USA.gov</a> ~" +
"    <a href='https://www.hhs.gov/vulnerability-disclosure-policy/index.html' class='footlink'>HHS Vulnerability Disclosure</a>" +
"  </div>" +
"</div>" +
"</div>" +
"<br />" +
"</body>" +
"</html>";

// print it out
document.write(footer);
