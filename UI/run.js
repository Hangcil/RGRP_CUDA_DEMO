function getProgressBar(testID) {
    switch (testID) {
        case 1:
            return progressBar1
        case 2:
            return progressBar2
        case 3:
            return progressBar3
        case 4:
            return progressBar4
        default:
            return
    }
}

function onGoClicked(testID) {
    app.inCalculating = true
    getProgressBar(testID).Material.accent = "#ff3b30"
    getProgressBar(testID).visible = true
    buttonGO.enabled = false
    buttonResult.enabled = false
    launcher.startSimulation(cflSlider.value, Number(endTInput.text), Number(gammaInput.text), alphaSlider.value, parseInt(cellNumXInput.text), parseInt(cellNumYInput.text), app.testID)
}

function onProgressChanged(testID, progress) {
    getProgressBar(testID).value = progress
}

function onSimulationCompleted() {
    app.inCalculating = false
    getProgressBar(testID).Material.accent = "#ff9500"
}

function onPlotCompleted(filePath) {
    buttonResult.enabled = true
    buttonGO.enabled = true
    app.inCalculating = false
    app.plotPath = filePath
    getProgressBar(testID).Material.accent = "#06943d"
}