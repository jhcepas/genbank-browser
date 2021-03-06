/* 
CT_HOST variable will be defined in the auto-generated config js file. 
YOu can overwrite it here:

var CT_HOST = 'http://localhost';
*/

var url = $.url();
var qregion = url.param('r');
var genomeViewer = null;

function EmblAdapter(args) {

    _.extend(this, Backbone.Events);

    _.extend(this, args);

    this.on(this.handlers);

    this.cache = {};
}

EmblAdapter.prototype = {

    getData: function (args) {
        var _this = this;

        args.webServiceCallCount = 0;


        var region = args.region;
        if (region.start > 300000000 || region.end < 1) {
            return;
        }
        region.start = (region.start < 1) ? 1 : region.start;
        region.end = (region.end > 300000000) ? 300000000 : region.end;


        var params = {};
        _.extend(params, this.params);
        _.extend(params, args.params);
        var dataType = args.dataType;
        if (_.isUndefined(dataType)) {
            console.log("dataType must be provided!!!");
        }
        var chunkSize;
        /********/


        if (dataType == 'histogram') {
            var histogramId = dataType + '_' + params.interval;
            if (_.isUndefined(this.cache[histogramId])) {
                this.cache[histogramId] = new FeatureChunkCache({chunkSize: params.interval});
            }
            chunkSize = this.cache[histogramId].chunkSize;
            // Extend region to be adjusted with the chunks
            var adjustedRegions = this.cache[histogramId].getAdjustedRegions(region);

            var url = CT_HOST + '/' + this.resource;

            if (adjustedRegions.length > 0) {
                params['region'] = adjustedRegions.toString();
                url = Utils.addQueryParamtersToUrl(params, url);
                $.ajax({
                    type: 'GET',
                    url: url,
                    success: function (data) {
                        _this._histogramSuccess(data, dataType, histogramId, args);
                    }
                });

            } else {
                var chunksByRegion = this.cache[histogramId].getCachedByRegion(region);
                var chunksCached = this.cache[histogramId].getByRegions(chunksByRegion.cached);
                this.trigger('data:ready', {items: chunksCached, dataType: dataType, chunkSize: chunkSize, sender: this});
            }


        } else {
            //Create one FeatureChunkCache by datatype
            if (_.isUndefined(this.cache[dataType])) {
                this.cache[dataType] = new FeatureChunkCache(this.cacheConfig);
            }
            chunkSize = this.cache[dataType].chunkSize;
            var chunksByRegion = this.cache[dataType].getCachedByRegion(region);

            if (chunksByRegion.notCached.length > 0) {
                var queryRegionStrings = _.map(chunksByRegion.notCached, function (region) {
                    return new Region(region).toString();
                });

                //limit queries
                var n = 50;
                var lists = _.groupBy(queryRegionStrings, function (a, b) {
                    return Math.floor(b / n);
                });
                var queriesList = _.toArray(lists); //Added this to convert the returned object to an array.

                for (var i = 0; i < queriesList.length; i++) {

                    /** CT custom **/
                    var url = CT_HOST + '/' + this.resource;
                    var queryParams = {
                        region: queriesList[i].toString().toLowerCase()
                    };
                    _.extend(queryParams, this.params);
                    url = Utils.addQueryParamtersToUrl(queryParams, url);
                    $.ajax({
                        type: 'GET',
                        url: url,
                        success: function (data) {
                            _this._success(data, dataType, args);
                        }
                    });
                }
            }
            if (chunksByRegion.cached.length > 0) {
                var chunksCached = this.cache[dataType].getByRegions(chunksByRegion.cached);
                this.trigger('data:ready', {items: chunksCached, dataType: dataType, chunkSize: chunkSize, sender: this});
            }
        }
        if (args.webServiceCallCount === 0) {
            args.done();
        }
    },

    _success: function (data, dataType, args) {
        var timeId = this.resource + " save " + Utils.randomString(4);
        /** time log **/

        var chunkSize = this.cache[dataType].chunkSize;

        var chunks = [];
        for (var i = 0; i < data.response.length; i++) {
            var queryResult = data.response[i];
            var region = new Region(queryResult.id);
            var features = queryResult.result;
            var chunk = this.cache[dataType].putByRegion(region, features);
            chunks.push(chunk);
        }

        if (chunks.length > 0) {
            this.trigger('data:ready', {items: chunks, dataType: dataType, chunkSize: chunkSize, sender: this});
        }
        if (args.webServiceCallCount === 0) {
            args.done();
        }


    },
    _histogramSuccess: function (data, dataType, histogramId, args) {

        var chunkSize = this.cache[histogramId].chunkSize;

        var chunks = [];
        for (var i = 0; i < data.response.length; i++) {
            var queryResult = data.response[i];
            for (var j = 0; j < queryResult.result.length; j++) {
                var interval = queryResult.result[j];
                var region = new Region(queryResult.id);
                region.load(interval);
                chunks.push(this.cache[histogramId].putByRegion(region, interval));
            }
        }
//        var chunksByRegion = this.cache[histogramId].getB(region);

        this.trigger('data:ready', {items: chunks, dataType: dataType, chunkSize: chunkSize, sender: this});
        if (args.webServiceCallCount === 0) {
            args.done();
        }
    }
};

function CtSequenceAdapter(args) {

    _.extend(this, Backbone.Events);

    this.id = Utils.genId("CtSequenceAdapter");

    //set default args
    this.host;
    this.gzip = true;

    //set instantiation args, must be last
    _.extend(this, args);

    this.sequence = {};
    this.start = {};
    this.end = {};

    this.on(this.handlers);
}

CtSequenceAdapter.prototype.clearData = function () {
    this.sequence = {};
    this.start = {};
    this.end = {};
};

CtSequenceAdapter.prototype.getData = function (args) {
    var _this = this;

    this.sender = args.sender;
    var region = args.region;
    var chromosome = region.chromosome;

    region.start = (region.start < 1) ? 1 : region.start;
    region.end = (region.end > 300000000) ? 300000000 : region.end;

    //clean when the new position is too far from current
    if (region.start < this.start[chromosome] - 5000 || region.end > this.end[chromosome] + 5000) {
        this.clearData();
    }
    console.log("seq ars,", args);
    var params = {};
    _.extend(params, this.params);


    var queryString = this._getSequenceQuery(region);

    if (queryString != "") {

        /** CT custom **/
        var url = CT_HOST + '/' + this.resource;
        var queryParams = {
            region: queryString.toLowerCase()
        };
        url = Utils.addQueryParamtersToUrl(queryParams, url);
        console.log(url);
        $.ajax({
            type: 'GET',
            url: url,
            success: function (data) {
                _this._processSequenceQuery(data, true);
            }
        });
        /**/

    } else {
        if (this.sender != "move") {
            // this.trigger('data:ready',{
            //     items:{
            //         sequence:this.sequence[chromosome],
            //         start:this.start[chromosome],
            //         end:this.end[chromosome]
            //     },
            //     params:params
            // });
            this.trigger('data:ready', {
                items: {
                    sequence: this.sequence[chromosome],
                    start: this.start[chromosome],
                    end: this.end[chromosome]
                },
                params: params,
                sender: this
            });
        }
    }

};

CtSequenceAdapter.prototype._getSequenceQuery = function (region) {
    var _this = this;
    var chromosome = region.chromosome;

    var s, e, query, querys = [];
    if (_this.start[chromosome] == null && _this.end[chromosome] == null) {
        //args.start -= 100;
        //args.end += 100;
        _this.start[chromosome] = region.start;
        _this.end[chromosome] = region.end;
        s = region.start;
        e = region.end;
        query = chromosome + ":" + s + "-" + e;
        querys.push(query);
    } else {
        if (region.start < _this.start[chromosome]) {
            s = region.start;
            e = _this.start[chromosome] - 1;
            e = (e < 1) ? region.end = 1 : e;
            _this.start[chromosome] = s;
            query = region.chromosome + ":" + s + "-" + e;
            querys.push(query);
        }
        if (region.end > _this.end[chromosome]) {
            e = region.end;
            s = _this.end[chromosome] + 1;
            _this.end[chromosome] = e;
            query = region.chromosome + ":" + s + "-" + e;
            querys.push(query);
        }
    }
    return querys.toString();
};

CtSequenceAdapter.prototype._processSequenceQuery = function (data, throwNotify) {
    var _this = this;
    var params = data.params;


    for (var i = 0; i < data.response.length; i++) {
        var queryResponse = data.response[i];
        var splitDots = queryResponse.id.split(":");
        var splitDash = splitDots[1].split("-");
        var queryStart = parseInt(splitDash[0]);
        var queryEnd = parseInt(splitDash[1]);

        var queryId = queryResponse.id;
        var seqResponse = queryResponse.result;
        var chromosome = seqResponse.chromosome;

        if (this.sequence[chromosome] == null) {
            this.sequence[chromosome] = seqResponse.sequence;
        } else {
            if (queryStart == this.start[chromosome]) {
                this.sequence[chromosome] = seqResponse.sequence + this.sequence[chromosome];
            } else {
                this.sequence[chromosome] = this.sequence[chromosome] + seqResponse.sequence;
            }
        }

        if (this.sender == "move" && throwNotify == true) {
            this.trigger('data:ready', {
                items: {
                    sequence: seqResponse.sequence,
                    start: queryStart,
                    end: queryEnd
                },
                params: params,
                sender: this
            });
        }
    }

    if (this.sender != "move" && throwNotify == true) {
        this.trigger('data:ready', {
            items: {
                sequence: this.sequence[chromosome],
                start: this.start[chromosome],
                end: this.end[chromosome]
            },
            params: params,
            sender: this
        });
    }
};

//Used by bam to get the mutations
CtSequenceAdapter.prototype.getNucleotidByPosition = function (args) {
    var _this = this;
    if (args.start > 0 && args.end > 0) {
        var queryString = this._getSequenceQuery(args);

        var chromosome = args.chromosome;

        if (queryString != "") {

            /** CT custom **/
            var url = CT_HOST + '/' + this.resource;
            var queryParams = {
                region: queryString.toLowerCase()
            };
            url = Utils.addQueryParamtersToUrl(queryParams, url);
            $.ajax({
                type: 'GET',
                url: url,
                async: false
            });
            _this._processSequenceQuery(data);
            /**/

        }
        if (this.sequence[chromosome] != null) {
            var referenceSubStr = this.sequence[chromosome].substr((args.start - this.start[chromosome]), 1);
            return referenceSubStr;
        } else {
            console.log("SequenceRender: this.sequence[chromosome] is undefined");
            return "";
        }
    }
};

var run = function () {
    /* region and species configuration */
    if (qregion) {
        var coords = qregion.split(",");

        var region = new Region({
            chromosome: coords[0],
            start: parseInt(coords[1]),
            end: parseInt(coords[2]),
        });
    } else {
        var region = new Region({
            chromosome: default_chr, 
            start: default_start,
            end: default_end,
        });
    }
    var availableSpecies = {
        "text": "Species",
        "items": [
            {
                "text": "Chaetomium",
                "items": [
                    {
                        "text": default_species,
                        "assembly": default_assembly,
                        "region": {
                            "chromosome": default_chr,
                            "start": default_start,
                            "end": default_end,
                        },
                        "chromosomes": scaffolds,
                        "url": url,
                    }
                ]
            }
        ]
    }
    var species = availableSpecies.items[0].items[0];

    /*** Quick search config ***/
    var quickSearchResultFn = function (queryStr) {
        var results;
        $.ajax({
            type: 'GET',
            url: CT_HOST + '/search/' + queryStr,
            dataType: 'json',
            async: false,
            success: function (data) {
                results = data.response.matches;
            }
        });
        return results;
    };
    /*****/


    genomeViewer = new GenomeViewer({
        target: 'application',
        width: $('#application').width(),
        region: region,
        availableSpecies: availableSpecies,
        species: species,
        sidePanel: false,
        resizable: true,
        drawKaryotypePanel: false,
        drawChromosomePanel: false,
        karyotypePanelConfig: {
            collapsed: false,
            collapsible: true
        },
        chromosomePanelConfig: {
            collapsed: false,
            collapsible: true
        },

        /*** Quick search config ***/
        quickSearchResultFn: quickSearchResultFn,
        quickSearchDisplayKey: '1',
        handlers: {
            'quickSearch:select': function (event) {
            }
        },
        navigationBarConfig: {
            componentsConfig: {
//                restoreDefaultRegionButton:false,
//                regionHistoryButton:false,
//                speciesButton:false,
//                chromosomesButton:false,
                karyotypeButton: false,
                chromosomeButton: false,
                regionButton: false,
//                zoomControl:false,
//                windowSizeControl:false,
//                positionControl:false,
//                moveControl:false,
//                autoheightButton:false,
//                compactButton:false,
                searchControl: false
            }
        },
        handlers: {
            'quickSearch:select': function (event) {
            }
        },
        version:'',
        chromosomeList: scaffolds_data
    });

    tracks = [];

    this.gccont = new FeatureTrack({
        title: 'GC content',
        /***/
        minHistogramRegionSize: 0,//nts
        maxLabelRegionSize: 10000000,//nts
        visibleRegionSize: 30000000000,//nts
        /***/
        histogramMaxFreqValue: 1,
        height: 50,
        renderer: new GeneRenderer(),

        dataAdapter: new EmblAdapter({
            host: CT_HOST,
            resource: "gc_content",
            params: {
            },
            species: genomeViewer.species,
            cacheConfig: {
                chunkSize: default_chunksize
            }
        })
    });
    tracks.push(this.gccont);

    this.sequence = new SequenceTrack({
        //title: 'Sequence',
        height: 23,
        visibleRegionSize: 200,//nts

        renderer: new SequenceRenderer(),

        dataAdapter: new CtSequenceAdapter({
            host: CT_HOST,
            subCategory: "region",
            resource: "sequence",
            species: genomeViewer.species,
            featureCache: {
                gzip: true,
                chunkSize: default_chunksize
            }
        })
    });
    tracks.push(this.sequence);

    this.gene = new GeneTrack({
        title: 'Genes (CDS, tRNA, rRNA, misc_RNA)',
        maxLabelRegionSize: 2000000000,
        minTranscriptRegionSize: 100000,
        histogramMaxFreqValue: 1,
        height: 400,

        renderer: new GeneRenderer(FEATURE_CONFIG.tRNA),

        dataAdapter: new EmblAdapter({
            host: CT_HOST,
            resource: "gene",
            params: {
                biotypes: 'CDS,tRNA,rRNA,misc_RNA,repeat_region,misc_feature'
            },
            species: genomeViewer.species,
            cacheConfig: {

                chunkSize: default_chunksize
            }
        })
    });
    this.gene.renderer.on({
        'feature:click': function (event) {
            show_gene_info(event);
        },
    });
    tracks.push(this.gene);

    var renderer = new FeatureRenderer(FEATURE_TYPES.gene);
    renderer.on({
        'feature:click': function (event) {
            show_gene_info(event);
        }
    });


    this.repeats = new FeatureTrack({
        title: 'Repeat Regions',
        maxLabelRegionSize: 999999999999,
        minTranscriptRegionSize: 10,
        histogramMaxFreqValue: 1,
        height: 400,

        renderer: renderer,

        dataAdapter: new EmblAdapter({
            host: CT_HOST,
            resource: "repeat_region",
            params: {
            },

            species: genomeViewer.species,
            cacheConfig: {
                chunkSize: default_chunksize
            }
        })
    });
    this.repeats.renderer.on({
        'feature:click': function (event) {
            show_gene_info(event);
        },
    });
    tracks.push(this.repeats);


    var geneOver = new FeatureTrack({
        //title: 'Region overview',
        maxLabelRegionSize: 999999999999,
        minTranscriptRegionSize: 1,
        //histogramZoom: 0,
        histogramMaxFreqValue: 1,
        //labelZoom: 40,
        height: 150,
        titleVisibility: 'hidden',
        featureTypes: FEATURE_TYPES,

        renderer: renderer,

        dataAdapter: new EmblAdapter({
            host: CT_HOST,
            resource: "gene",
            params: {
                exclude: 'transcripts',
            },
            species: genomeViewer.species,
            cacheConfig: {
                chunkSize: default_chunksize
            }
        })
    });

    genomeViewer.addOverviewTrack(geneOver);
    genomeViewer.addTrack(tracks);

    genomeViewer.draw();


    $('#positionControl').css({width: '330px' });
    $('#autoheightButton').trigger('click');
};

jQuery(document).ready(function () {
    jQuery("#seqid_search").select2(
        {placeholder: "gene ID or description term",
            multiple: false,
            width: '500px',
            allowClear: true,
            dropdownAutoWidth: true,
            containerCss: {"border": "0px"},
            containerCssClass: "search_field",
            closeOnSelect: false,
            openOnEnter: true,
            minimumInputLength: 3,
            formatSearching: 'Searching...',
            formatNoMatches: 'Sorry, no matches found',
            ajax: {
                type: 'post',
                url: CT_HOST + '/search/',
                data: function (term, page, context) {
                    return {'seqid': term};
                },
                dataType: 'json',
                results: function (data, page) {
                    return data;
                },
            },

            initSelection: function (element, callback) {
                $.ajax({
                    type: 'post',
                    url: CT_HOST + '/search/',
                    data: {'seqid': element.val()},
                    dataType: 'json',
                    success: function (data, page) {
                        console.log(data.results[0]);
                        callback(data.results[0]);
                    },
                })
            },

            formatResult: function (e, container, query) {
                var desc = e.desc;
                if (e.desc.length > 10000) {
                    desc = e.desc.substring(0, 120) + "...";
                }
                var uniprot = '';
                if (e.uniprot.length > 0) {
                    var uniprot = " (" + e.uniprot + ")";
                }

                html = '<span class="s_geneid">' + e.id + '</span>' +
                    '<span class="s_uniprot">' + uniprot + '</span>' +
                    '<span class="s_biotype">, ' + e.biotype + '. </span>' +
                    '<span class="s_desc"> - ' + desc + ' -</span>' +
                    '<span class="s_region"> [' + e.reg + ']</span>';
                var regex = new RegExp('(' + query.term + ')', 'gi');
                var hi_html = html.replace(regex, '<span class="hi">$1</span>');

                return hi_html;


            },
            formatSelection: function (e) {
                html = '<span class="s_geneid">' + e.id + '</span>' +
                    '<span class="s_uniprot"> (' + e.uniprot + ')</span>' +
                    '<span class="s_biotype"> - ' + e.biotype + ' -</span>' +
                    '<span class="s_region"> [' + e.reg + ']</span>';
                return html;
            },
            escapeMarkup: function (m) {
                return m;
            },
        });

    $('#seqid_search').on('change', function (event) {
        update_search();
    });
    $('#seqid_search').select2('focus');

});

function sleep(millis, callback, args) {
    setTimeout(function () {
        callback(args);
    }, millis);
}

function update_search() {
    var locus = $('#seqid_search').select2('data');
    
    if (locus != undefined) {
        console.log('setting region', locus.chr, parseInt(locus.start) - 2000, parseInt(locus.end) + 2000);
        genomeViewer.setRegion({chromosome: locus.chr, 
            start: parseInt(locus.start) - 2000,
            end: parseInt(locus.end) + 2000});
    }
}

function go_to_gene(genename) {
    $("#ph_genes").popover('hide');
    $("#expr_genes").popover('hide');
    console.log('got to gene', genename)
    $("#seqid_search").select2('val', genename);
    sleep(250, update_search);
}

function show_gene_info(gene) {
    console.log(gene.feature.id);

    $('#speciesButton').popover({
        html: true,
        placement: "bottom",
        trigger: 'click',
        title: '<span id="gene_ext_title" style="margin-right:20px;">' + gene.feature.id + '</span><button type="button" id="close" class="close" onclick="$(&quot;#speciesButton&quot;).popover(&quot;hide&quot;);">&times;</button>',
        content: '<div id="gene_ext_info">Loading... </div>',
        container: "body",
    });

    $('#speciesButton').popover('show');
    $.ajax({
        type: 'POST',
        url: CT_HOST + '/seqid_data/',
        data: {seqid: gene.feature.fid},
        success: function (data) {
            var h = '';
            for (var key in {'type': 0, 'location': 0, 'qualifiers': 0, 'strand': 0,
                'pep_img': 0, 'length:': 0, 'NT Seq': 0}) {
                if (key in data) {
                    if (key == 'pep_img') {
                        h += '<b><a target="_blank" href="http://wlab.ethz.ch/protter/">Protter</a> peptide image:</b><br> <a target="_blank" href="' + data[key] + '"><img height="300" src="' + data[key] + '"></a>'
                    }
                    else if (key == 'qualifiers') {
                        var trans = '';
                        for (var qkey in data[key]) {
                            if (qkey == 'translation') {
                                trans = '<b>' + qkey + ': </b><pre>' + data[key][qkey] + '</pre>';
                            } else {
                                h += '<b>' + qkey + ': </b>' + data[key][qkey] + '<br>';
                            }
                        }
                        h += trans;
                    }
                    else if (key == 'NT Seq') {
                        h += '<b>nt seq:</b><pre>' + data[key] + '</pre>';
                    }

                    else {
                        h += '<b>' + key + ': </b>' + data[key] + '<br>';
                    }
                }
            }
            $('#gene_ext_info').html(h);
        }
    });

}

