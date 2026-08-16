document.addEventListener("DOMContentLoaded", function(event) {
    const banner = document.getElementById('banner');
    if (typeof banner !== 'undefined' && banner !== null) {
        const header = document.querySelector('header.d-print-none');
        if (typeof header !== 'undefined' && header !== null) {
            header.setAttribute('id', 'has-banner');
        }
    }

    const homeHeader = document.getElementById('headers-content-container');
    if (typeof homeHeader !== 'undefined' && homeHeader !== null) {
        const headerH1 = document.querySelector('.header-get-function');
        if (typeof headerH1 !== 'undefined' && headerH1 !== null) {
            if(headerH1.classList.contains('home-header')){
                homeHeader.classList.add('home-header-container');
            }
        }
    }
});

$(function () {
    if(window.isAdBlockActive){
        let modalAdBlockerWarning = $('.modal.adblocker-warning');
        if(modalAdBlockerWarning.length > 0){
            $(modalAdBlockerWarning).modal('show');
        }
    }


  // Remove any previously injected trigger
  var existing = document.querySelector('.search-mobile-trigger');
  if (existing) existing.remove();

  var searchContainer = document.querySelector('#search_content_container');
  var breadcrumbsContainer = document.querySelector('#breadcrumbs_content_container');
  if (!searchContainer || !breadcrumbsContainer) return;

  searchContainer.classList.remove('search-expanded');

  // Create trigger button
  var triggerDiv = document.createElement('div');
  triggerDiv.className = 'search-mobile-trigger';

  var btn = document.createElement('button');
  btn.setAttribute('title', 'Start site search');
  btn.className = 'btn btn-sm btn-success';
  btn.innerHTML = '<span class="co-search"></span>';
  triggerDiv.appendChild(btn);

  breadcrumbsContainer.insertAdjacentElement('afterend', triggerDiv);

  btn.addEventListener('click', function (e) {
    e.stopPropagation();
    var isExpanded = searchContainer.classList.contains('search-expanded');
    if (!isExpanded) {
      // Show search row 2, hide trigger button
      searchContainer.classList.add('search-expanded');
      triggerDiv.classList.add('search-trigger-hidden');
      var input = document.querySelector('#search_query_solr');
      if (input) setTimeout(function () { input.focus(); }, 50);
    } else {
      // Hide search row 2, restore trigger button
      searchContainer.classList.remove('search-expanded');
      triggerDiv.classList.remove('search-trigger-hidden');
    }
  });

  // Close on outside click — register only once
  if (!window._searchMobileDocListener) {
    window._searchMobileDocListener = function (e) {
      if (window.innerWidth >= 400) return;
      var sc = document.querySelector('#search_content_container');
      var td = document.querySelector('.search-mobile-trigger');
      if (!sc || !td) return;
      if (!sc.contains(e.target) && !td.contains(e.target)) {
        sc.classList.remove('search-expanded');
        td.classList.remove('search-trigger-hidden');
      }
    };
    document.addEventListener('click', window._searchMobileDocListener);
  }
});