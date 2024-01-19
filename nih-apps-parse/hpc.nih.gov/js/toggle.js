// This is the gross way of toggling content
$(document).ready(function(){

// When the page is loaded, all the toggle-content elements are hidden
//    $("*").load(function() {
//        $(".toggle-tag").each(function(i){
//            $(this).siblings(".toggle-content").toggle();
          //$(this).siblings(".toggle-content").css("display: none");
// Keep the old format, for the most part
//            $(this).css( "color", "#0000d5" );
//            $(this).css( "font-weight", "bold" );
//            $(this).css( "font-size", "16px" );
//            $(this).css( "padding-left", "15px" );
//            $(this).css( "padding-top", "8px" );
//            $(this).css("text-decoration","none") 
//            $(this).mouseover(function() { 
//                $(this).css("text-decoration","underline") 
//                $(this).css("cursor","hand") 
//            });
//            $(this).mouseout(function() { 
//                $(this).css("text-decoration","none") 
//                $(this).css("cursor","pointer") 
//            });
//        });
//    });

// When an element of class toggle-icon is clicked, toggle the icon class and all elements of class toggle-content within the toggleable shell
    $(".toggle-tag").click(function(){
        $("i",this).toggleClass("fa-caret-down fa-caret-right");
        $(this).siblings(".toggle-content").toggle();
    });


// Globally toggle all the toggle-icons and toggle-contents
    $(".toggle-all").click(function(){

// Toggle the icon for the global toggle switch
        $("i",this).toggleClass("fa-caret-right fa-caret-down");

// Hide everything
        if ($("i",this).is(".fa-caret-right")) {
            $(".toggle-content").each(function(i){
                $(this).hide();
            });
            $("i.fa-caret-down").each(function(i){
                $(this).removeClass("fa-caret-down");
                $(this).addClass("fa-caret-right");
            });
        }
// Show everything
        else {
            $(".toggle-content").each(function(i){
                $(this).show();
            });
            $("i.fa-caret-right").each(function(i){
                $(this).removeClass("fa-caret-right");
                $(this).addClass("fa-caret-down");
            });
        }
    });
    
    // in addition to the toggle above, here is a simple
    // tab view. Styling in main.css
    $('ul.tabs li').click(function(){
        var tab_id = $(this).attr('data-tab');

        $('ul.tabs li').removeClass('current');
        $('.tab-content').removeClass('current');

        $(this).addClass('current');
        $("#"+tab_id).addClass('current');
    });
    // tab view END
});
