FEATURE_CONFIG = {
    gene: {
        filters: [
            {
                name: "biotype",
                text: "Biotype",
                values: ["3prime_overlapping_ncrna", "ambiguous_orf", "antisense", "disrupted_domain", "IG_C_gene", "IG_D_gene", "IG_J_gene", "IG_V_gene", "lincRNA", "miRNA", "misc_RNA", "Mt_rRNA", "Mt_tRNA", "ncrna_host", "nonsense_mediated_decay", "non_coding", "non_stop_decay", "polymorphic_pseudogene", "processed_pseudogene", "processed_transcript", "protein_coding", "pseudogene", "retained_intron", "retrotransposed", "rRNA", "sense_intronic", "sense_overlapping", "snoRNA", "snRNA", "transcribed_processed_pseudogene", "transcribed_unprocessed_pseudogene", "unitary_pseudogene", "unprocessed_pseudogene"],
                selection: "multi"
            }
        ]
        //options:[
        //]
    },

};
FEATURE_OPTIONS = {
    gene: [
        {
            name: "biotype",
            text: "Biotype",
            values: ["3prime_overlapping_ncrna", "ambiguous_orf", "antisense", "disrupted_domain", "IG_C_gene", "IG_D_gene", "IG_J_gene", "IG_V_gene", "lincRNA", "miRNA", "misc_RNA", "Mt_rRNA", "Mt_tRNA", "ncrna_host", "nonsense_mediated_decay", "non_coding", "non_stop_decay", "polymorphic_pseudogene", "processed_pseudogene", "processed_transcript", "protein_coding", "pseudogene", "retained_intron", "retrotransposed", "rRNA", "sense_intronic", "sense_overlapping", "snoRNA", "snRNA", "transcribed_processed_pseudogene", "transcribed_unprocessed_pseudogene", "unitary_pseudogene", "unprocessed_pseudogene"],
            selection: "multi"
        }
    ],

};

GENE_BIOTYPE_COLORS = {
    "gene": "seagreen",
    "CDS": "seagreen",
    "tRNA": "red",
    "rRNA": "#8b668b",
    "misc_RNA": "YellowGreem",
    "hypothetical protein": "yellow",
    "repeat_region": "steelblue",
    "repeat_region2": "orange",
    "exon": "yellow",
    null: "black"
};


SEQUENCE_COLORS = {A: "#009900", C: "#0000FF", G: "#857A00", T: "#aa0000", N: "#555555"};


FEATURE_TYPES = {

    //methods
    formatTitle: function (str) {
        var s = '';
        if(str){
            str.replace(/_/gi, " ");
            s = s.charAt(0).toUpperCase() + s.slice(1);
        }
        return s;
    },
    getTipCommons: function (f) {
        var strand = (f.strand != null) ? f.strand : "NA";
        return 'start-end:&nbsp;<span class="emph">' + f.start + '-' + f.end + '</span><br>' +
            'strand:&nbsp;<span class="emph">' + strand + '</span><br>' +
            'length:&nbsp;<span class="info">' + (f.end - f.start + 1).toString().replace(/(\d)(?=(\d\d\d)+(?!\d))/g, "$1,") + '</span><br>';
    },

    //items
    sequence: {
        color: SEQUENCE_COLORS
    },
    undefined: {
        getLabel: function (f) {
            var str = "";
            str += f.chromosome + ":" + f.start + "-" + f.end;
            return str;
        },
        getTipTitle: function (f) {
            return " ";
        },
        getTipText: function (f) {
            return " ";
        },
        getColor: function (f) {
            return "grey";
        },
//		infoWidgetId: "id",
        height: 10,
        histogramColor:"blue"
    },
    gene: {
        label: function (f, zoom) {
            var name = (f.name != null) ? f.name : f.id;
            var str = "";
            str += (f.strand < 0 || f.strand == '-') ? "<" : "";
            str += " " + name + " ";
            str += (f.strand > 0 || f.strand == '+') ? ">" : "";
            //if (f.biotype != null && f.biotype != '' && f.biotype != 'gene' && zoom > 25) {
            //    str += " [" + f.biotype + "]";
            //}
            return str;
        },
        tooltipTitle: function (f) {
            return FEATURE_TYPES.formatTitle('Gene') +
                ' - <span class="">' + f.id + ' ('+f.uniprot+')</span>';
        },


        tooltipText: function (f) {
            var color = GENE_BIOTYPE_COLORS[f.biotype];
            var html = '<div style="width:500px;"> id:&nbsp;<span class="ssel">' + f.id + '</span><br>' +
                'biotype:&nbsp;<span class="emph" style="color:' + color + ';">' + f.biotype + '</span><br>' +
                FEATURE_TYPES.getTipCommons(f);

            for (var key in f.annotations) {
                html += "<b>" + key+ "</b>:&nbsp;" + f.annotations[key] + "<br>";}
            for (var key in f.expression) {
                html += "<b>" + key+ "</b>:&nbsp;" + f.annotations[key] + "<br>";}

            if (f.pep_img && f.pep_img == 1){
                html+= '<img width="280px" src="pep_img/png/' +f.id
                    + '.png"><br>';
            }
            return html+"</div>";
        },
        color: function (f) {
            return GENE_BIOTYPE_COLORS[f.biotype];
        },

        infoWidgetId: "id",
        height: function(f){
            if (f.biotype == "gene"){
                return 5;}
            else {
                return 8;}
        },
        histogramColor: "lightblue"
    },

    transcript: {
        label: function (f, zoom) {
            var name = f.id;
            var str = "";
            str += (f.strand < 0) ? "<" : "";
            str += " " + name + " ";
            str += (f.strand > 0) ? ">" : "";
            if (f.biotype != null && f.biotype != '') {
                str += " [" + f.biotype + "]";
            }
            return str;
        },
        tooltipTitle: function (f) {
            return FEATURE_TYPES.formatTitle('Transcript') +
                ' - <span class="ok">' + f.id + ' ('+f.uniprot+')</span>';
        },
        tooltipText: function (f) {
            var color = GENE_BIOTYPE_COLORS[f.biotype];
            return    'id:&nbsp;<span class="ssel">' + f.id + '</span><br>' +
                'biotype:&nbsp;<span class="emph" style="color:' + color + ';">' + f.biotype + '</span><br>' +
                'product:&nbsp;<span class="emph">' + f.product + '</span><br>' +
                'description:&nbsp;<span class="emph">' + f.description + '</span><br>' +
                FEATURE_TYPES.getTipCommons(f);
        },
        color: function (f) {
            return "black";
        },
        infoWidgetId: "id",
        height: 1,
        histogramColor: "lightblue"
    },
    exon: {
        label: function (f) {
            var name = (f.name != null) ? f.name : f.id;
            return name;
        },
        tooltipTitle: function (f) {
            var name = (f.name != null) ? f.name : f.id;
            if (name == null) {
                name = ''
            }
            return FEATURE_TYPES.formatTitle('Exon') + ' - <span class="ok">' + name + '</span>';
        },
        tooltipText: function (e, t) {
            var ename = (e.name != null) ? e.name : e.id;
            var tname = (t.name != null) ? t.name : t.id;
            var color = "black";
            return    'transcript name:&nbsp;<span class="ssel">' + t.name + '</span><br>' +
                'transcript Ensembl&nbsp;ID:&nbsp;<span class="ssel">' + t.id + '</span><br>' +
                'transcript biotype:&nbsp;<span class="emph" style="color:' + color + ';">' + t.biotype + '</span><br>' +
                'product:&nbsp;<span class="emph">' + t.product + '</span><br>' +
                'transcript description:&nbsp;<span class="emph">' + t.description + '</span><br>' +
                'transcript start-end:&nbsp;<span class="emph">' + t.start + '-' + t.end + '</span><br>' +
                'exon start-end:&nbsp;<span class="emph">' + e.start + '-' + e.end + '</span><br>' +
                'strand:&nbsp;<span class="emph">' + t.strand + '</span><br>' +
                'length:&nbsp;<span class="info">' + (e.end - e.start + 1).toString().replace(/(\d)(?=(\d\d\d)+(?!\d))/g, "$1,") + '</span><br>';
        },
        color: function (f) {
            return 'orange';
        },
        infoWidgetId: "id",
        height: 4,
        histogramColor: "lightblue"
    },

};
