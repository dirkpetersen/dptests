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
document.write("Last modified: " + date + " " + lmonth + " " + fyear + " ");
