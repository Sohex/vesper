/*
 *
 */
var CoPublisher = CoPublisher || {};
CoPublisher.Helper = CoPublisher.Helper || {};
CoPublisher.Helper.Color = CoPublisher.Helper.Color || (function($) {
        'use strict';
//privates

        var _color;

        var parseColorCode =function parseColor(color) {
            var cleanedColor = $.trim(color);
            var firstChar = cleanedColor.charAt(0);
            var rgb = {red: 0, green: 0, blue: 0};

            if(firstChar == '#')
            {
                var colorOfNumber = cleanedColor.substr(1,cleanedColor.length - 1);
                var red, green, blue = 0;

                if(colorOfNumber.length < 6)
                {
                    red = colorOfNumber.slice(0,1);
                    red.concat(red);

                    green = colorOfNumber.slice(1,2);
                    green.concat(green);

                    blue = colorOfNumber.slice(2,3);
                    blue.concat(blue);

                }
                else
                {
                    red = colorOfNumber.slice(0,2);
                    green = colorOfNumber.slice(2,4);
                    blue = colorOfNumber.slice(4,6);
                }
                rgb.red = parseInt(red,16);
                rgb.green = parseInt(green,16);
                rgb.blue = parseInt(blue,16);
            }
            if(firstChar.toLowerCase() == 'r')
            {
                cleanedColor = cleanedColor.toLowerCase();
                cleanedColor = cleanedColor.replace(/rgb/gi, "");
                cleanedColor = cleanedColor.replace(/\(/gi, "");
                cleanedColor = cleanedColor.replace(/\)/gi, "");

                var colorArray = cleanedColor.split(",");

                if(colorArray[0] != undefined)
                    rgb.red = parseInt($.trim(colorArray[0]),10);
                if(colorArray[1] != undefined)
                    rgb.green = parseInt($.trim(colorArray[1]),10);
                if(colorArray[2] != undefined)
                    rgb.blue = parseInt($.trim(colorArray[2]),10);
                console.log("it's a rgb color ");
                console.log(rgb);
            }

            if(firstChar.toLowerCase() == 'h')
                throw "Only hex and rgb color support";

            return rgb;
        }
        var lightColorCode = function lightenColor(percentage) {

            percentage = (percentage/100);

            var rgb = parseColorCode(_color);
            rgb.red = Math.round(rgb.red * (1 - percentage));
            rgb.green = Math.round(rgb.green * (1 - percentage));
            rgb.blue = Math.round(rgb.blue * (1 - percentage));

            return rgb;
        };

//public interface
        return {
            init: function(htmlColorCode) {
                _color = htmlColorCode;
            },
            lighten: function(percentage){
                var rgb = lightColorCode(percentage);
                return "rgb("+rgb.red+","+rgb.green+","+rgb.blue+")";
            }
        };

    }) (jQuery);
CoPublisher.JournalMetrics = CoPublisher.JournalMetrics || {};
CoPublisher.JournalMetrics.Model = CoPublisher.JournalMetrics.Model || (function() {

        'use strict';

        //privates
        var _monthAsString = ["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct", "Nov", "Dec"];
        var _categories;
        var _series;
        var _windowWidth;
        var _max;

        var calculateMaxNumberOfBars = function calculateMaxNumberOfBars() {
            //console.log("sWidth "+window.innerWidth);
            var width = _windowWidth || window.innerWidth;
            return width < 1024 ? 8 : 9;
        };

        var amountOfElementsToStockUp = function calculateAmountOfElementsToStockUp(categoriesLength){

            if(categoriesLength < _max)
            {
                var currentMonthCount = categoriesLength;
                return _max - currentMonthCount;
            }

            return 0;
        };

        var calculateDateToStart = function calculateDateToStart(categoriesLength){

            var now = new Date(Date.now());
            var thisMonthAsString = now.getFullYear()+ ' ' + _monthAsString[now.getMonth()];
            var lastDate = categoriesLength > 0 ? _categories[categoriesLength-1] : thisMonthAsString;

            var splitedDate = lastDate.split(' ',2);
            return new Date(splitedDate[0]+" 01, "+splitedDate[1]);

        };

        var generateData = function generateSeriesAndCategories(last,maxRuns) {
            var newElement;

            for(var index = 1; index <= maxRuns ; index ++) {

                last.setMonth(last.getMonth() + 1);
                //console.log(index);
                newElement = _monthAsString[last.getMonth()]+" "+last.getFullYear();
                _categories.push(newElement);
                _series[0]['data'].push(0);
                _series[1]['data'].push(0);
                _series[2]['data'].push(0);
            }
        };

        //public interface
        return {
            init: function(categories,series,windowWidth) {
                _categories = categories;
                _series = series;
                _windowWidth = windowWidth;
                _max = calculateMaxNumberOfBars();
                var categoriesLength = _categories.length;

                var numberOfElementsToStockUp = amountOfElementsToStockUp(categoriesLength);
                var last = calculateDateToStart(categoriesLength);
                generateData(last,numberOfElementsToStockUp);

            },
            getSeries: function() {
                return _series;
            },
            getCategories: function() {
                return _categories;
            },
            getXAxisMin: function() {
                return _categories.length - _max;
            },
            getXAxisMax: function() {
                return _categories.length - 1;
            }
        };

    })();

CoPublisher.JournalMetrics.View = CoPublisher.JournalMetrics.View || (
        function ($) {

            'use strict';

            var _chartSelector;
            var _viewModel;
            var _default = {
                mainColor: '#606060',
                borderColor: '#DEDEDE',
                local: 'en-GB'
            };
            var _settings;

            var renderHighChart = function renderJournalMetricBarChart() {

                var highChartInit = function() {
                    //console.log('load highchart');
                    var chart = this, counter = 0;

                    $.each(chart.legend.allItems, function(i, item){
                        var $check = $(item.checkbox),
                            left = parseFloat($check.css('left')),
                            label = item.legendItem,
                            parentID = $check.parent().attr('id'),
                            checkBoxID = parentID+'--category-selector-'+counter;

                        $check.css({
                            left: (left - item.checkboxOffset + 18 ) + 'px'
                        }).attr('checked','checked').attr('id',checkBoxID);

                        $('<label/>')
                            .attr('for',checkBoxID)
                            .css('position',$check.css('position'))
                            .css('width','16')
                            .css('height',$check.css('height'))
                            .css('left',$check.css('left'))
                            .css('top',$check.css('top'))
                            .addClass('checked')
                            .insertBefore( $check )
                            .on( "click", function() {

                                var labelID = $(this).attr('for'),
                                                res = labelID.split("--"),
                                                parentId = res[0],
                                                elem = res[1].split('-'),
                                                childNum = parseInt(elem[elem.length-1], 10),
                                                rectSelector = '#'+parentId+' .highcharts-legend-item rect'; 

                                $(rectSelector).each(function(index, item){

                                    if(index == childNum)
                                    {
                                        console.log('trigger click on rest '+index);
                                        $(item).trigger('click');

                                        return false;
                                    }

                                });

                        } );

                        counter++;
                    });
                };

                var toLocalString = function (number) {

                    if(typeof number === 'number' && typeof number.toLocaleString === "function" )
                        return number.toLocaleString(_settings.local);

                    return number;
                };

                var legendItemClickedActions = function(event) {
                    var clickedCheckBox = event.currentTarget.checkbox;
                    var checkedCurrent = !$(clickedCheckBox).prop('checked');

                    $(clickedCheckBox).prop('checked', checkedCurrent);

                    var $fakeCheckBox = $(event.currentTarget.checkbox).prev();

                    var isChecked = $fakeCheckBox.hasClass('checked');

                    if(isChecked)
                        $fakeCheckBox.removeClass('checked');
                    else
                        $fakeCheckBox.addClass('checked');
                };

                $(_chartSelector).highcharts({
                    credits: false,
                    scrollbar: {
                        enabled: true,
                        barBackgroundColor: _settings.mainColor,
                        barBorderRadius: 3,
                        barBorderWidth: 0,
                        barBorderColor: _settings.borderColor,
                        buttonBackgroundColor: '#FFF',
                        buttonBorderWidth: 0,
                        buttonBorderRadius: 0,
                        buttonArrowColor: _settings.mainColor,
                        rifleColor: _settings.mainColor,
                        trackBackgroundColor: 'none',
                        trackBorderWidth: 0,
                        trackBorderRadius: 0,
                        trackBorderColor: '#fff'
                    },
                    chart: {
                        type: 'column',
                        marginBottom: 65,
                        events: {
                            load: highChartInit
                        }
                    },
                    colors:    ['#69CD4B', '#DC5046', '#398CDA'],
                    title: {
                        text: ''
                    },
                    xAxis: {
                        categories: _viewModel.getCategories(),
                        min: _viewModel.getXAxisMin(),
                        max: _viewModel.getXAxisMax()
                    },
                    yAxis: {
                        min: 0,
                        title: {
                            text: ''
                        },
                        stackLabels: {
                            enabled: true,
                            style: {
                                fontWeight: 'bold',
                                color: _settings.mainColor
                            },
                            formatter: function() {
                                if(this.total > 0)
                                    return toLocalString(this.total);

                                return null;
                            }
                        }
                    },
                    legend: {
                        align: 'center',
                        x: 0,
                        verticalAlign: 'bottom',
                        y: 10,
                        floating: false,
                        backgroundColor: 'white',
                        borderColor: '#CCC',
                        borderWidth: 1,
                        shadow: false,
                        itemStyle: {
                            fontWeight: 'normal',
                            color: _settings.mainColor
                        },
                        itemHiddenStyle: {
                            color: _settings.mainColor
                        }
                    },
                    tooltip: {
                        borderRadius: 0,
                        formatter: function() {
                            return '<b>'+ this.x +'</b><br/>'+
                                this.series.name +': '+ toLocalString(this.y) +'<br/>'+
                                'Total: '+ toLocalString(this.point.stackTotal);
                        }
                    },
                    plotOptions: {
                        line: {
                            animation: false
                        },
                        series: {
                            showCheckbox:true,
                            selected:true,
                            animation: false,
                            dataLabels:{
                                enabled:true,
                                color: 'white',
                                formatter:function(){

                                        return toLocalString(this.y);
                                },
                                style: {
                                    textShadow: '1px 1px 1px #000'
                                },
                                overflow: false,
                                inside: true
                            },
                            events: {
                                legendItemClick: legendItemClickedActions
                            }
                        },
                        column: {
                            animation: true,
                            stacking: 'normal'
                        }


                    },
                    series: _viewModel.getSeries()
                });

                /*
                 Resize Arrow Dimension
                 d = M 8 4 L 8 10 5 7 => M 9 2 L 9 12 3 7
                 d = M 6 4 L 6 10 9 7 => M 5 2 L 5 12 11 7

                 I'm to lazy to parse
                 */
                $('.highcharts-scrollbar > g > path').each(function(index, value){

                    var dimension = $(value).attr('d');
                    switch(dimension)
                    {
                        case 'M 8 4 L 8 10 5 7':
                            $(value).attr('d','M 9 2 L 9 12 3 7');
                            break;
                        case 'M 6 4 L 6 10 9 7':
                            $(value).attr('d','M 5 2 L 5 12 11 7');
                            break;
                    }
                });

            };


//public interface
            return {
                init: function(chartContainer,model,config) {
                    _chartSelector = chartContainer;
                    _viewModel = model;
                    _settings = $.extend({}, _default, config);
                },
                render: function(){
                    renderHighChart();
                }
            };
        }) (jQuery);//import