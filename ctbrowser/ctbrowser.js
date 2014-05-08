CELLBASE_HOST = ""
CT_HOST = "http://ctbrowser.embl.de";

var scaffolds =   ["scf7180000011816",
                   "scf7180000011820",
                   "scf7180000011822",
                   "scf7180000011806",
                   "scf7180000011818",
                   "scf7180000011814",
                   "scf7180000011821",
                   "scf7180000011815",
                   "scf7180000011819",
                   "scf7180000011812",
                   "scf7180000011823",
                   "scf7180000011801",
                   "scf7180000011817",
                   "rrn",
                   "scf7180000011813",
                   "scf7180000011802",
                   "scf7180000011824",
                   "scf7180000011808",
                   "scf7180000011809",
                   "scf7180000011807",
                   "scf7180000011810",
                   "mito"]; 

var scaffolds_data = [
    {
        end: 8133,
        name: "rrn",
        cytobands: [ ],
        start: 1,
        numberGenes: 0,
        isCircular: 0,
        size: 8133
    },
    {
        end: 1603,
        name: "scf7180000011824",
        cytobands: [ ],
        start: 1,
        numberGenes: 0,
        isCircular: 0,
        size: 1603
    },
    {
        end: 1013,
        name: "scf7180000011808",
        cytobands: [ ],
        start: 1,
        numberGenes: 0,
        isCircular: 0,
        size: 1013
    },
    {
        end: 1150,
        name: "scf7180000011809",
        cytobands: [ ],
        start: 1,
        numberGenes: 0,
        isCircular: 0,
        size: 1150
    },
    {
        end: 3461039,
        name: "scf7180000011820",
        cytobands: [ ],
        start: 1,
        numberGenes: 0,
        isCircular: 0,
        size: 3461039
    },
    {
        end: 1965831,
        name: "scf7180000011821",
        cytobands: [ ],
        start: 1,
        numberGenes: 0,
        isCircular: 0,
        size: 1965831
    },
    {
        end: 3123527,
        name: "scf7180000011822",
        cytobands: [ ],
        start: 1,
        numberGenes: 0,
        isCircular: 0,
        size: 3123527
    },
    {
        end: 745067,
        name: "scf7180000011823",
        cytobands: [ ],
        start: 1,
        numberGenes: 0,
        isCircular: 0,
        size: 745067
    },
    {
        end: 13765,
        name: "scf7180000011802",
        cytobands: [ ],
        start: 1,
        numberGenes: 0,
        isCircular: 0,
        size: 13765
    },
    {
        end: 192410,
        name: "scf7180000011801",
        cytobands: [ ],
        start: 1,
        numberGenes: 0,
        isCircular: 0,
        size: 192410
    },
    {
        end: 2849635,
        name: "scf7180000011806",
        cytobands: [ ],
        start: 1,
        numberGenes: 0,
        isCircular: 0,
        size: 2849635
    },
    {
        end: 11746,
        name: "scf7180000011807",
        cytobands: [ ],
        start: 1,
        numberGenes: 0,
        isCircular: 0,
        size: 11746
    },
    {
        end: 127206,
        name: "mito",
        cytobands: [ ],
        start: 1,
        numberGenes: 0,
        isCircular: 0,
        size: 127206
    },
    {
        end: 1047547,
        name: "scf7180000011819",
        cytobands: [ ],
        start: 1,
        numberGenes: 0,
        isCircular: 0,
        size: 1047547
    },
    {
        end: 2847644,
        name: "scf7180000011818",
        cytobands: [ ],
        start: 1,
        numberGenes: 0,
        isCircular: 0,
        size: 2847644
    },
    {
        end: 1853679,
        name: "scf7180000011815",
        cytobands: [ ],
        start: 1,
        numberGenes: 0,
        isCircular: 0,
        size: 1853679
    },
    {
        end: 2150193,
        name: "scf7180000011814",
        cytobands: [ ],
        start: 1,
        numberGenes: 0,
        isCircular: 0,
        size: 2150193
    },
    {
        end: 153753,
        name: "scf7180000011817",
        cytobands: [ ],
        start: 1,
        numberGenes: 0,
        isCircular: 0,
        size: 153753
    },
    {
        end: 6909506,
        name: "scf7180000011816",
        cytobands: [ ],
        start: 1,
        numberGenes: 0,
        isCircular: 0,
        size: 6909506
    },
    {
        end: 1124,
        name: "scf7180000011810",
        cytobands: [ ],
        start: 1,
        numberGenes: 0,
        isCircular: 0,
        size: 1124
    },
    {
        end: 10547,
        name: "scf7180000011813",
        cytobands: [ ],
        start: 1,
        numberGenes: 0,
        isCircular: 0,
        size: 10547
    },
    {
        end: 973894,
        name: "scf7180000011812",
        cytobands: [ ],
        start: 1,
        numberGenes: 0,
        isCircular: 0,
        size: 973894
    }
]

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
                "assembly": "EMBL v1.0",
                "region": {
                    "chromosome":"scf7180000011816",
                    "start": 379031,
                    "end": 404027,
                },
                "chromosomes": scaffolds,
                
                "url": "ct.bork.embl.de"
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
            console.log(event, "click in transcript");
        }
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


