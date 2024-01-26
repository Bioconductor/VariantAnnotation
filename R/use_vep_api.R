#' helper function to construct inputs for VEP REST API
#' @import httr
#' @import jsonlite
#' @param chr character(1)
#' @param pos numeric(1)
#' @param id character(1)
#' @param ref character(1)
#' @param alt character(1)
#' @note Produces a string used as an example in VEP API documentation.
variant_body = function(chr, pos, id, ref, alt) {
  sprintf("%s  %d  %s %s %s . . .", chr, pos, id, ref, alt)
}
#' elementary vep/homo_sapiens/region call to ensembl VEP REST API
#' @importFrom curl has_internet
#' @param chr character(1) ensembl chromosome identifier (e.g., "7")
#' @param pos numeric(1) 1-based chromosome position
#' @param id character(1) arbitrary identifier
#' @param ref character(1) reference allele
#' @param alt character(1) alternative allele
#' @note This function prepares a POST to rest.ensembl.org/vep/homo_sapiens/region endpoint.
#' @return Instance of 'response' defined in httr package.
#' @examples
#' chk = post_Hs_region("7", 155800001, "chk", "A", "T")
#' chk
#' res = jsonlite::fromJSON(jsonlite::toJSON(httr::content(chk)))
#' dim(chk)
#' @export
post_Hs_region = function(chr, pos, id, ref, alt) {
 stopifnot(curl::has_internet())
 server <- "https://rest.ensembl.org"
 ext <- "/vep/homo_sapiens/region"
 tj = jsonlite::toJSON(list(variants=variant_body(chr, pos, id, ref, alt)))
 ans = httr::POST(paste(server, ext, sep = ""), 
   httr::content_type("application/json"), 
   httr::accept("application/json"), 
   body = as.character(tj)
   )
 httr::stop_for_status(ans)
 ans
}

#' Use the VEP region API on variant information in a VCF object as defined in VariantAnnotation.
#' @param vcfobj instance of VCF class; note the difference between the CollapsedVCF and
#' ExpandedVCF instances.
#' @param snv_only logical(1) if TRUE filter the VCF to information about single nucleotide addresses
#' @param chk_max logical(1) requests to ensembl VEP API are limited to 200 positions; if
#' TRUE and the request involves more than 200 positions, an error is thrown by this function.
#' @return instance of 'response' from httr package
#' @examples
#' fl <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
#' r22 = readVcf(fl)
#' dr = which(width(rowRanges(r22))!=1)
#' r22s = r22[-dr]
#' res = vep_by_region(r22[1:100], snv_only=FALSE, chk_max=FALSE)
#' ans = jsonlite::fromJSON(jsonlite::toJSON(httr::content(res)))
#' @export
vep_by_region = function(vcfobj, snv_only=TRUE, chk_max=TRUE) {
 stopifnot(inherits(vcfobj, "VCF"))
 if (snv_only) {
  dr = which(width(rowRanges(vcfobj))!=1)
  if (length(dr)>0) vcfobj = vcfobj[-dr]
 }
 if (chk_max) {
  if (nrow(vcfobj) > 200) stop("VEP API is limited to 200 positions")
  }
 rr = rowRanges(vcfobj)
 post_Hs_region( chr = as.character(seqnames(rr)),
                 pos = start(rr),
                 id = names(rr),
                 ref = as.character(rr$REF),
                 alt = as.character(unlist(rr$ALT)) )
}


