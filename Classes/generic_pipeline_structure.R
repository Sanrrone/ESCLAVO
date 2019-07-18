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



