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
    snp: [
        {
            name: "consequence_type",
            text: "Consequence Type",
            values: ["2KB_upstream_variant", "5KB_upstream_variant", "500B_downstream_variant", "5KB_downstream_variant", "3_prime_UTR_variant", "5_prime_UTR_variant", "coding_sequence_variant", "complex_change_in_transcript", "frameshift_variant", "incomplete_terminal_codon_variant", "inframe_codon_gain", "inframe_codon_loss", "initiator_codon_change", "non_synonymous_codon", "intergenic_variant", "intron_variant", "mature_miRNA_variant", "nc_transcript_variant", "splice_acceptor_variant", "splice_donor_variant", "splice_region_variant", "stop_gained", "stop_lost", "stop_retained_variant", "synonymous_codon"],
            selection: "multi"
        }
    ],
    bam: [
        {
            name: "view",
            text: "View",
            values: ["view_as_pairs", "show_soft-clipped_bases"],
            selection: "multi"
        }
    ]
};

GENE_BIOTYPE_COLORS = {
    "gene": "seagreen",
    "CDS": "seagreen",
    "tRNA": "red",
    "rRNA": "#8b668b",
    "misc_RNA": "YellowGreem",
    "hypothetical protein": "yellow",
    "repeat_region": "steelblue",
    null: "black"
};

// 
//    "antisense": "SteelBlue",
//    "disrupted_domain": "YellowGreen",
//    "IG_C_gene": "#FF7F50",
//    "IG_D_gene": "#FF7F50",
//    "IG_J_gene": "#FF7F50",
//    "IG_V_gene": "#FF7F50",
//    "lincRNA": "#8b668b",
//    "miRNA": "#8b668b",
//    "misc_RNA": "#8b668b",
//    "Mt_rRNA": "#8b668b",
//    "Mt_tRNA": "#8b668b",
//    "ncrna_host": "Fuchsia",
//    "nonsense_mediated_decay": "seagreen",
//    "non_coding": "orangered",
//    "non_stop_decay": "aqua",
//    "polymorphic_pseudogene": "#666666",
//    "processed_pseudogene": "#666666",
//    "processed_transcript": "#0000ff",
//    "protein_coding": "#a00000",
//    "pseudogene": "#666666",
//    "retained_intron": "goldenrod",
//    "retrotransposed": "lightsalmon",
//    "rRNA": "indianred",
//    "sense_intronic": "#20B2AA",
//    "sense_overlapping": "#20B2AA",
//    "snoRNA": "#8b668b",
//    "snRNA": "#8b668b",
//    "transcribed_processed_pseudogene": "#666666",
//    "transcribed_unprocessed_pseudogene": "#666666",
//    "unitary_pseudogene": "#666666",
//    "unprocessed_pseudogene": "#666666",
//    "": "orangered",
//    "other": "#000000"
//};


SNP_BIOTYPE_COLORS = {
    "2KB_upstream_variant": "#a2b5cd",
    "5KB_upstream_variant": "#a2b5cd",
    "500B_downstream_variant": "#a2b5cd",
    "5KB_downstream_variant": "#a2b5cd",
    "3_prime_UTR_variant": "#7ac5cd",
    "5_prime_UTR_variant": "#7ac5cd",
    "coding_sequence_variant": "#458b00",
    "complex_change_in_transcript": "#00fa9a",
    "frameshift_variant": "#ff69b4",
    "incomplete_terminal_codon_variant": "#ff00ff",
    "inframe_codon_gain": "#ffd700",
    "inframe_codon_loss": "#ffd700",
    "initiator_codon_change": "#ffd700",
    "non_synonymous_codon": "#ffd700",
    "intergenic_variant": "#636363",
    "intron_variant": "#02599c",
    "mature_miRNA_variant": "#458b00",
    "nc_transcript_variant": "#32cd32",
    "splice_acceptor_variant": "#ff7f50",
    "splice_donor_variant": "#ff7f50",
    "splice_region_variant": "#ff7f50",
    "stop_gained": "#ff0000",
    "stop_lost": "#ff0000",
    "stop_retained_variant": "#76ee00",
    "synonymous_codon": "#76ee00",
    "other": "#000000"
};


SEQUENCE_COLORS = {A: "#009900", C: "#0000FF", G: "#857A00", T: "#aa0000", N: "#555555"};

SAM_FLAGS = [
    ["read paired", 0x1],
    ["read mapped in proper pair", 0x2],
    ["read unmapped", 0x4],
    ["mate unmapped", 0x8],
    ["read reverse strand", 0x10],
    ["mate reverse strand", 0x20],
    ["first in pair", 0x40],
    ["second in pair", 0x80],
    ["not primary alignment", 0x100],
    ["read fails platform/vendor quality checks", 0x200],
    ["read is PCR or optical duplicate", 0x400]
];


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
                html+= '<img width="280px" src="pep_img/svg/' +f.id
                    + '.svg" onerror="this.src=\'pep_img/png/'+f.id+'.png\';"><br>'; 
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
            return GENE_BIOTYPE_COLORS[f.product];
        },
        infoWidgetId: "id",
        height: 4,
        histogramColor: "lightblue"
    },

};

