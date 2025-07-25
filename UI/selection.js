function setConfig(text, Nx, Ny, T, gamma, CFL, alpha, imageSource) {
    dispTestNameRectText.text = text;
    cellNumXInput.text = Nx
    cellNumYInput.text = Ny
    endTInput.text = T
    gammaInput.text = gamma
    cflSlider.value = CFL
    alphaSlider.value = alpha
    preview.source = imageSource
}

function updateConfig() {
    buttonResult.enabled = false
    plotPath = ""
    progressBar1.visible = false
    progressBar2.visible = false
    progressBar3.visible = false
    progressBar4.visible = false
    switch (app.testID) {
        case 3:
            setConfig(qsTr("Mach 800 Jet Problem"), 600, 200, 0.002, 1.4, 0.80, 1.90, "../images/jet800.png")
            break;
        case 1:
            setConfig(qsTr("Interaction of 4 Shocks"), 200, 200, 0.80, 1.4, 0.50, 2, "../images/4shocks.png")
            break;
        case 2:
            setConfig(qsTr("Interaction of 4 Contact Discontinuities"), 200, 200, 0.80, 1.4, 0.50, 2, "../images/4cds.png")
            break;
        case 4:
            setConfig(qsTr("Double Mach Reflection Problem"), 800, 200, 0.2, 1.4, 0.80, 2, "../images/2Mach.png")
            break;
        default:
            break;
    }

}