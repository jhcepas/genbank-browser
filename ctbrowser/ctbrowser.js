CT_HOST = "http://ctbrowser.embl.de";

var url = $.url();
var qregion = url.param('r');
var genomeViewer = null;
var run = function() {
    /* region and species configuration */
    if (qregion){
        var coords = qregion.split(",");
        console.log(coords.length);
        
        var region = new Region({
            chromosome: coords[0],
            start: parseInt(coords[1]),
            end: parseInt(coords[2]),
        });
    } else {
        var region = new Region({
            chromosome: "scf7180000011816",
            start: 379031,
            end: 404027,
        });
    }
    console.log(region);
    var availableSpecies= {
        "text": "Species",
        "items": [{
            "text": "Chaetomium",
            "items": [{
                "text": "Chaetomium thermophilum",
                "assembly": "EMBL v2.0",
                "region": {
                    "chromosome":"scf7180000011816",
                    "start": 379031,
                    "end": 404027,
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
			type:'GET',
			url:CT_HOST+'/search/'+queryStr, //tRNA
			dataType: 'json',
			async:false,
			success:function(data){
				results = data.response.matches;
			}
		});
		return results;
	};
	/*****/
    
    genomeViewer = new GenomeViewer({
        targetId: 'application',
        region: region,
        availableSpecies: availableSpecies,
        species: species,
        sidePanel: true,
        autoRender: true,
        border: true,
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
		handlers:{
			'quickSearch:select': function (event) {
				console.log(event.item);
			}
		},
        
        chromosomeList: scaffolds_data,
        
		/*****/
    }); //the div must exist
    
    genomeViewer.draw();
    
    tracks = [];
    
    this.gccont = new FeatureTrack({
        targetId: null,
        id: 3,
        title: 'GC content',
        /***/
        minHistogramRegionSize: 0,//nts
		maxLabelRegionSize: 10000000,//nts
		visibleRegionSize:30000000000,//nts
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
                chunkSize: 500000
            }
        })
    });
    
    tracks.push(this.gccont);
    
    
    
    this.sequence = new SequenceTrack({
        targetId: null,
        id: 1,
        //title: 'Sequence',
        height: 20,
        visibleRegionSize:300,//nts
        
        renderer: new SequenceRenderer(),
        
        dataAdapter: new CtSequenceAdapter({
            host: CT_HOST,
            subCategory: "region",
            resource: "sequence",
            species: genomeViewer.species,
            featureCache: {
                gzip: true,
                chunkSize: 1000
            }
        })
    });
    
    tracks.push(this.sequence);
    
    this.gene = new GeneTrack({
        targetId: null,
        id: 2,
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
                biotypes: 'CDS,tRNA,rRNA,misc_RNA'
            },
            
            species: genomeViewer.species,
            cacheConfig: {
                chunkSize: 500000
            }
        })
    });
    this.gene.renderer.on({
        'feature:click': function(event) {
            console.log(event, "click in transcript"); },
        'feature:click': function(event) {
            console.log(event, "click in transcript"); }



    });
    tracks.push(this.gene);
    
    var renderer = new FeatureRenderer('gene');
    renderer.on({
        'feature:click': function(event) {
            console.log(event, "click in gene");
        }
    });
    var geneOver = new FeatureTrack({
        targetId: null,
        id: 4,
        //title: 'Region overview',
        maxLabelRegionSize: 200000000,
		minTranscriptRegionSize: 100000,
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
                chunkSize: 500000
            }
        })
    });
    genomeViewer.addOverviewTrack(geneOver);
    
    genomeViewer.addTrack(tracks);
    
    
	/*custom search*/
	var customSearchButton = $('#customSearchButton')[0];
	var customSearchField = $('#customSearchField')[0];
	var getCustomSearchResults = function(queryStr){
		var html = '<div id="customSearchDiv">';
		var results = [];
		$.ajax({
			type:'GET',
			url:CT_HOST+'/search/'+queryStr, //tRNA
			dataType: 'json',
			async:false,
			success:function(data){
				results = data.response.matches;
			}
		});
        html += "<table>"
        
		for(var i = 0;i<results.length;i++){
			var result = results[i];
            html += "<tr style='border-bottom: 1px solid grey;'><td>";
            html += '<span class="customSearchItem" style="over" id="' + result[1] + '" region="'+result[0] +'" >' + result[1] + '</span>' ;
            html += '</td><td>';
            html += '<span class="customSearchItemType" style="over" id="' + result[1] + '" region="'+result[0] +'" >[' + result[2] + ']</span>' ;
            html += '</td><td>';
            html += '<span class="customSearchItemDesc" style="over" id="' + result[1] + '" region="'+result[0] +'" >' + result[3] + '</span>' ;
            html += '</td><tr>';
		}
		html += '</table></div>';
		return html;
	};
    
	$(customSearchButton).click(function(){
		var queryStr = $(customSearchField).val();
		if($('.popover').length > 0){
			$(customSearchButton).popover('hide');
		}else{
			if(queryStr !== ''){
				$(customSearchButton).popover('destroy');
				var html = getCustomSearchResults(queryStr);
				$(customSearchButton).popover({
					title:'Search results:',
					html:true,
					content:html,
					placement:'bottom',
					trigger:'manual', 
                    
				});
				$(customSearchButton).popover('show');
				$('#customSearchDiv').click(function(event){
				    $(customSearchButton).popover('destroy');
					var regionStr = $(event.target).attr('region');
					console.log(event.target)
					console.log(this);
					var region = new Region(regionStr);
					genomeViewer.setRegionDirect(region);
                    
				});
                
			}
		}
        
	});
	
};

jQuery(document).ready(function(){
    jQuery("#seqid_search").select2(
        {placeholder: "gene name",
         multiple:false,
         width:'50%',
         allowClear:true,
         dropdownAutoWidth:true, 
         containerCss:{"border":"0px"},
         containerCssClass:"search_field",
         closeOnSelect: false,
         openOnEnter: true,
         minimumInputLength: 3,
         formatSearching: 'Searching...',
         formatNoMatches: 'Sorry, no matches found',
         
         ajax: { 
             type: 'post',             
             url: CT_HOST+'/search/',
             data: function (term, page, context){ return {'seqid': term}; }, 
             dataType: 'json',
             results: function (data, page) { 
                 console.log(data);
                 return data;
             },
         },
         formatResult: function(e, container, query) {
             var regex = new RegExp( '(' + query.term + ')', 'gi' );
             var high = e.text.replace(regex, '<b class="hi2">$1</b>');
             return '<span>'+e.id+'</span>'+'<span> '+e.biotype+'</span>'+'<span>'+e.desc+'</span>'+'<span>'+e.reg+'</span>'


         },
         formatSelection: function(e) {
                var sp_string = '<span class="selectable_sp" onClick="javascript:show_target_species(this, \''+e.id+'\');">'+e.species+'</span>';
             return e.oid + " from " + sp_string;
         },
         escapeMarkup: function (m) { return m; },
        });
    
});