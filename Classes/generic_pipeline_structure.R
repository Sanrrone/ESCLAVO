setClass("pipeline", representation=list(
  id="character",
  name="character",
  version="numeric",
  image="character",
  steps="list",
  stepsOptions="list"
  
  )
)
setGeneric("getAnalysisID", function(ob) {standardGeneric("getAnalysisID")})
setMethod("getAnalysisID", "pipeline", function(ob) {return(ob@id)})

setGeneric("getAnalysisName", function(ob) {standardGeneric("getAnalysisName")})
setMethod("getAnalysisName", "pipeline", function(ob) {return(ob@name)})

setGeneric("getAnalysisVersion", function(ob) {standardGeneric("getAnalysisVersion")})
setMethod("getAnalysisVersion", "pipeline", function(ob) {return(ob@version)})

setGeneric("getAnalysisSteps", function(ob) {standardGeneric("getAnalysisSteps")})
setMethod("getAnalysisSteps", "pipeline", function(ob) {return(ob@steps)})

setGeneric("getAnalysisStepsNames", function(ob) {standardGeneric("getAnalysisStepsNames")})
setMethod("getAnalysisStepsNames", "pipeline", function(ob) {return(names(ob@steps))})

setGeneric("getStep", function(ob,step) {standardGeneric("getStep")})
setMethod("getStep", "pipeline", function(ob,step) {return(ob@steps[[step]])})

setGeneric("setNewStepFolder", function(ob,step,newfolder) {standardGeneric("setNewStepFolder")})
setMethod("setNewStepFolder", "pipeline", function(ob,step,newfolder) {
  ob@steps[[step]]$folder<-newfolder
  return(ob)
})


setGeneric("getAnalysisOptions", function(ob) {standardGeneric("getAnalysisOptions")})
setMethod("getAnalysisOptions", "pipeline", function(ob) {
          accordionItem(
              id = ob@id,
              title = ob@name,
              color = ob@stepsOptions[["color"]],
              collapsed = ob@stepsOptions[["collapsed"]],
              "This is some text!"
            )
              
       }
)

setGeneric("getAnalysisImages", function(ob) {standardGeneric("getAnalysisImages")})
setMethod("getAnalysisImages", "pipeline", function(ob) {
  userListItem(
    src = dataURI(file=ob@image, mime="image/png"),
    user_name = ob@name, 
    description = "Updated: August 2019"
  )
})



