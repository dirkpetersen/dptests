//var host = window.location.origin;
var host = 'https://hpc.nih.gov';
var path = window.location.pathname;

if ((path === '/') || (path === '/index.php')) {
  document.write("<!DOCTYPE html PUBLIC '-//W3C//DTD XHTML 1.0 Transitional//EN' 'http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd'>" +
"<html xmlns='http://www.w3.org/1999/xhtml'>" +
"<head>" +
"  <meta http-equiv='X-UA-Compatible' content='IE=edge'>" +
"  <meta http-equiv='Content-Type' content='text/html; charset=utf-8' />" +
"  <script async type='text/javascript' id='_fed_an_ua_tag' src='https://dap.digitalgov.gov/Universal-Federated-Analytics-Min.js?agency=HHS'></script>" +
"  <link href='/css/main.css'  rel='stylesheet' type='text/css' />" +
"  <link href='/css/dropdown.css' rel='stylesheet' type='text/css' />" +
"  <script src='/js/jquery-3.5.1.min.js' type='text/javascript'></script>" +
"  <script src='/js/dropdown.js?v=" + Math.floor(Date.now())  + "' type='text/javascript'></script>" +
"  <link href='/css/font-awesome.min.css' rel='stylesheet' type='text/css' />" +
"  <link href='/css/quicklinks-dropdown.css' rel='stylesheet' type='text/css' />" +
"  <script src='/js/toggle.js' type='text/javascript'></script>" +
"  <!--[if lte IE 7.]>" +
"    <script defer type='text/javascript' src='/js/pngfix.js' ></script>" +
"  <![endif]-->" +
"  <title>NIH HPC Systems</title>" +
"</head>" +
"<body>" +
"<br />" +
"<div id='skiplink'><a href='#content-area' id='skiplink'></a></div>" +
"  <script type='text/javascript' language='JavaScript' src='" + host + "/js/banner.js?v=" + Math.floor(Date.now()) + "'></script>" +
"<div class='container'>" +
"  <table class='head_area_top' cellpadding='0' cellspacing='0' border='0' width='100%'>" +
"    <tr>" +
"      <th scope='col'></th>" +
"      <td valign='bottom' scope='col'>" +
"        <div class='titleLocation'>" +
"          <a href='/'><img src='" + host + "/images/BIOWULF_LOGO' border=0 height=72 alt='Biowulf High Performance Computing at the NIH'></a>" +
"        </div>" +
"      </td>" +
"      <td align='right' valign='top' scope='col'>" +
"        <script type='text/javascript' language='JavaScript' src='" + host + "/js/searchbox.js'></script>" +
"      </td>" +
"    </tr>" +
"  </table>" +
"  <script type='text/javascript' language='JavaScript' src='" + host + "/js/menubar.js?v=" + Math.floor(Date.now()) + "'></script>" +
"  <div class='main'>" +
"    <a name='content-area'></a>");
}
else {
  document.write("<!DOCTYPE html PUBLIC '-//W3C//DTD XHTML 1.0 Transitional//EN' 'http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd'>" +
"<html xmlns='http://www.w3.org/1999/xhtml'>" +
"<head>" +
"  <meta http-equiv='X-UA-Compatible' content='IE=edge'>" +
"  <meta http-equiv='Content-Type' content='text/html; charset=utf-8' />" +
"<!-- force browsers to reload -->" +
"  <meta http-equiv='Cache-Control' content='no-cache, no-store, must-revalidate' />" +
"  <meta http-equiv='Pragma' content='no-cache' />" +
"  <meta http-equiv='Expires' content='0' />" +
"  <link href='/css/main.css' rel='stylesheet' type='text/css' />" +
"  <link href='/css/dropdown.css' rel='stylesheet' type='text/css' />" +
"  <script async type='text/javascript' id='_fed_an_ua_tag' src='https://dap.digitalgov.gov/Universal-Federated-Analytics-Min.js?agency=HHS'></script>" +
"  <script src='/js/jquery-3.5.1.min.js' type='text/javascript'></script>" +
"  <script src='/js/dropdown.js?v=" + Math.floor(Date.now()) + "' type='text/javascript'></script>" +
"  <link href='/css/font-awesome.min.css' rel='stylesheet' type='text/css' />" +
"  <link href='/css/quicklinks-dropdown.css' rel='stylesheet' type='text/css' />" +
"  <script src='/js/toggle.js' type='text/javascript'></script>" +
"  <!--[if lte IE 7.]>" +
"    <script defer type='text/javascript' src='/js/pngfix.js' ></script>" +
"  <![endif]-->" +
"  <title>TEMPLATE</title>" +
"</head>" +
"<body>" +
"  <br />" +
"  <div id='skiplink'><a href='#content-area' id='skiplink'></a></div>" +
"  <script type='text/javascript' language='JavaScript' src='" + host + "/js/banner.js?v=" + Math.floor(Date.now()) + "'></script>" +
"  <div class='container'>" +
"    <table class='head_area' cellpadding='0' cellspacing='0' border='0' width='100%'>" +
"      <tr><th scope='col'></th>" +
"      <td valign='bottom' scope='col'>" +
"        <div class='titleLocation'>" +
"          <a href='" + host + "'><img src='/images/BIOWULF_LOGO' border=0 height=72 alt='Biowulf High Performance Computing at the NIH'></a>" +
"        </div>" +
"      </td>" +
"      <td align='right' valign='top' scope='col'>" +
"        <script type='text/javascript' language='JavaScript' src='" + host + "/js/searchbox.js'></script>" +
"      </td>" +
"    </tr>" +
"  </table>" +
"  <script type='text/javascript' language='JavaScript' src='" + host + "/js/menubar.js?v=" + Math.floor(Date.now()) + "'></script>" +
"  <div class='main'>" +
"    <a name='content-area'></a>");
}
