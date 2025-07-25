import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import QtQuick.Controls.Basic
import QtQuick.Controls.Material
import "selection.js" as UpdateConfig
import "run.js" as KernelRelated
import "data.js" as DataProcess

ApplicationWindow {
    id: app
    width: 800
    height: 600
    flags: Qt.FramelessWindowHint | Qt.Window
    title: "RGRP_CUDA"
    color: "transparent"
    visible: true
    property int spacing: width / 40 //20
    property int testID: 1
    property bool inCalculating: false
    property string plotPath: ""
    property int dragX: 0
    property int dragY: 0
    property bool dragging: false
    Rectangle {
        id: mainWind
        anchors.fill: parent
        color: "#ececec"
        radius: 5

        Rectangle {
            id: title
            x: 2.0 * app.spacing
            y: 0.15 * app.spacing
            width: 1.5 * app.spacing
            height: app.spacing
            color: "#00FFFFFF"
            Text {
                x: 0
                y: 0
                text: qsTr("RGRP CUDA DEMO")
                font.pointSize: 10
                font.family: "Microsoft YaHei"
            }
        }

        Rectangle {
            id: logo
            x: 0.2 * app.spacing
            y: title.y * 1.1
            width: 1.5 * app.spacing
            height: app.spacing
            color: "#00FFFFFF"
            Image {
                anchors.fill: parent
                anchors.margins: 2
                //source: "images/logo.png"
                source: "../images/logo.png"
                fillMode: Image.PreserveAspectFit
                horizontalAlignment: Image.AlignHCenter
                verticalAlignment: Image.AlignVCenter
                smooth: true
            }
        }

        MouseArea {
            width: mainWind.width
            height: 5 * app.spacing
            onPressed: {
                app.dragX = mouseX;
                app.dragY = mouseY;
                app.dragging = true;
            }
            onReleased: app.dragging = false
            onPositionChanged: {
                if (app.dragging) {
                    app.x += mouseX - app.dragX;
                    app.y += mouseY - app.dragY;
                }
            }
        }

        Rectangle {
            id: buttonClose
            x: mainWind.width - 1.65 * app.spacing
            y: 0.15 * app.spacing
            width: 1.5 * app.spacing
            height: app.spacing
            color: "#ececec"
            radius: 5
            Text {
                x: app.spacing / 3.5
                y: -app.spacing / 1.8
                text: qsTr("×")
                font.pointSize: 20
            }
            MouseArea {
                anchors.fill: parent
                hoverEnabled: true
                onEntered: parent.color = "#9c9b9b"
                onExited: parent.color = "#ececec"
                onClicked: app.close()
            }
        }
        Rectangle {
            id: buttonMinimize
            x: mainWind.width - 3.3 * app.spacing
            y: 0.15 * app.spacing
            width: 1.5 * app.spacing
            height: app.spacing
            color: "#ececec"
            radius: 5
            Text {
                x: app.spacing / 2.4
                y: -app.spacing / 1.5
                text: qsTr("-")
                font.pointSize: 25
            }
            MouseArea {
                anchors.fill: parent
                hoverEnabled: true
                onEntered: parent.color = "#9c9b9b"
                onExited: parent.color = "#ececec"
                onClicked: app.showMinimized()
            }
        }

        Rectangle {
            id: shocks4Button
            x: parent.width / 40
            y: parent.height / 30 + 0.5 * app.spacing
            width: parent.width * 0.3
            height: parent.height / 7.5
            radius: app.spacing / 2
            Text {
                id: shocks4ButtonText
                anchors.centerIn: parent
                font.pointSize: 20
                text: qsTr("Preset 1")
            }
            ProgressBar {
                id: progressBar1
                anchors.left: parent.left
                anchors.right: parent.right
                anchors.bottom: parent.bottom
                anchors.leftMargin: 2 * app.spacing
                anchors.rightMargin: 2 * app.spacing
                anchors.bottomMargin: app.spacing / 2
                Material.accent: "#ff3b30"
                visible: false
            }
            MouseArea {
                anchors.fill: parent
                hoverEnabled: true
                ToolTip.visible: containsMouse
                ToolTip.text: qsTr("interaction of 4 planar shocks")
                ToolTip.delay: 100
                onEntered: {
                    if (!app.inCalculating) {
                        parent.color = "#c1c1c1";
                    }
                }
                onExited: {
                    if (!app.inCalculating) {
                        parent.color = "white";
                    }
                }
                onClicked: {
                    if (!app.inCalculating) {
                        app.testID = 1;
                        UpdateConfig.updateConfig();
                    }
                }
            }
        }

        Rectangle {
            id: cds4Button
            x: shocks4Button.x
            y: shocks4Button.y + parent.height / 6
            width: shocks4Button.width
            height: shocks4Button.height
            radius: shocks4Button.radius
            Text {
                anchors.centerIn: parent
                font.pointSize: shocks4ButtonText.font.pointSize
                text: qsTr("Preset 2")
            }
            ProgressBar {
                id: progressBar2
                anchors.left: parent.left
                anchors.right: parent.right
                anchors.bottom: parent.bottom
                anchors.leftMargin: 2 * app.spacing
                anchors.rightMargin: 2 * app.spacing
                anchors.bottomMargin: app.spacing / 2
                Material.accent: "#ff3b30"
                visible: false
            }
            MouseArea {
                anchors.fill: parent
                hoverEnabled: true
                ToolTip.visible: containsMouse
                ToolTip.text: qsTr("interaction of 4 contact discontinuities")
                ToolTip.delay: 100
                onEntered: {
                    if (!app.inCalculating) {
                        parent.color = "#c1c1c1";
                    }
                }
                onExited: {
                    if (!app.inCalculating) {
                        parent.color = "white";
                    }
                }
                onClicked: {
                    if (!app.inCalculating) {
                        app.testID = 2;
                        UpdateConfig.updateConfig();
                    }
                }
            }
        }

        Rectangle {
            id: jet800Button
            x: shocks4Button.x
            y: cds4Button.y + parent.height / 6
            width: shocks4Button.width
            height: shocks4Button.height
            radius: shocks4Button.radius
            Text {
                anchors.centerIn: parent
                font.pointSize: shocks4ButtonText.font.pointSize
                text: qsTr("Preset 3")
            }
            ProgressBar {
                id: progressBar3
                anchors.left: parent.left
                anchors.right: parent.right
                anchors.bottom: parent.bottom
                anchors.leftMargin: 2 * app.spacing
                anchors.rightMargin: 2 * app.spacing
                anchors.bottomMargin: app.spacing / 2
                Material.accent: "#ff3b30"
                visible: false
            }
            MouseArea {
                anchors.fill: parent
                hoverEnabled: true
                ToolTip.visible: containsMouse
                ToolTip.text: qsTr("Mach 800 jet")
                ToolTip.delay: 100
                onEntered: {
                    if (!app.inCalculating) {
                        parent.color = "#c1c1c1";
                    }
                }
                onExited: {
                    if (!app.inCalculating) {
                        parent.color = "white";
                    }
                }
                onClicked: {
                    if (!app.inCalculating) {
                        app.testID = 3;
                        UpdateConfig.updateConfig();
                    }
                }
            }
        }

        Rectangle {
            id: mach2Button
            x: shocks4Button.x
            y: jet800Button.y + parent.height / 6
            width: shocks4Button.width
            height: shocks4Button.height
            radius: shocks4Button.radius
            Text {
                anchors.centerIn: parent
                font.pointSize: shocks4ButtonText.font.pointSize
                text: qsTr("Preset 4")
            }
            ProgressBar {
                id: progressBar4
                anchors.left: parent.left
                anchors.right: parent.right
                anchors.bottom: parent.bottom
                anchors.leftMargin: 2 * app.spacing
                anchors.rightMargin: 2 * app.spacing
                anchors.bottomMargin: app.spacing / 2
                Material.accent: "#ff3b30"
                visible: false
            }
            MouseArea {
                anchors.fill: parent
                hoverEnabled: true
                ToolTip.visible: containsMouse
                ToolTip.text: qsTr("double Mach reflection")
                ToolTip.delay: 100
                onEntered: {
                    if (!app.inCalculating) {
                        parent.color = "#c1c1c1";
                    }
                }
                onExited: {
                    if (!app.inCalculating) {
                        parent.color = "white";
                    }
                }
                onClicked: {
                    if (!app.inCalculating) {
                        app.testID = 4;
                        UpdateConfig.updateConfig();
                    }
                }
            }
        }

        Rectangle {
            id: previewArea
            x: shocks4Button.x
            y: mach2Button.y + parent.height / 6
            width: shocks4Button.width
            height: shocks4Button.height + app.spacing * 4
            radius: shocks4Button.radius - 2
            color: "#cfdeee"
            MouseArea {
                anchors.fill: parent
                hoverEnabled: true
                ToolTip.visible: containsMouse
                ToolTip.text: qsTr("preview")
                ToolTip.delay: 100
            }
            Image {
                id: preview
                anchors.fill: parent
                anchors.margins: 5
                source: "../images/4shocks.png"
                fillMode: Image.PreserveAspectFit
                horizontalAlignment: Image.AlignHCenter
                verticalAlignment: Image.AlignVCenter
                smooth: true
            }
        }

        Rectangle {
            id: configRectangle
            x: shocks4Button.x + shocks4Button.width + app.spacing
            y: shocks4Button.y
            width: parent.width / 1.6
            height: mach2Button.y + mach2Button.height + app.spacing
            color: "white"
            radius: shocks4Button.radius
            Rectangle {
                id: dispTestNameRect
                x: app.spacing
                y: app.spacing / 2
                width: configRectangle.width - x
                height: configRectangle.height / 20
                color: "white"
                radius: shocks4Button.radius
                Text {
                    id: dispTestNameRectText
                    text: qsTr("Interaction of 4 Shocks")
                    font.pointSize: shocks4ButtonText.font.pointSize
                }
            }
            TextField {
                id: cellNumXInput
                x: 2.5 * app.spacing
                y: dispTestNameRect.y + 1.5 * app.spacing + dispTestNameRect.height
                width: configRectangle.width / 2.5 - app.spacing
                height: configRectangle.height / 5 - app.spacing
                validator: IntValidator {
                    bottom: 0
                    top: 10000
                }
                text: qsTr("200")
                font.pointSize: shocks4ButtonText.font.pointSize - 5
                font.family: "Consolas"
                placeholderTextColor: "#d3d3d3"
                placeholderText: qsTr("Cell Number X")
                color: "black"
                background: Rectangle {
                    color: "#f5f5f5"
                    radius: app.spacing / 3
                }
            }
            TextField {
                id: cellNumYInput
                x: cellNumXInput.width + 4 * app.spacing
                y: cellNumXInput.y
                width: configRectangle.width / 2.5 - app.spacing
                height: configRectangle.height / 5 - app.spacing
                validator: IntValidator {
                    bottom: 0
                    top: 10000
                }
                text: qsTr("200")
                font.family: cellNumXInput.font.family
                font.pointSize: shocks4ButtonText.font.pointSize - 5
                placeholderTextColor: "#d3d3d3"
                placeholderText: qsTr("Cell Number Y")
                color: "black"
                background: Rectangle {
                    color: "#f5f5f5"
                    radius: app.spacing / 3
                }
            }
            TextField {
                id: endTInput
                x: cellNumXInput.x
                y: cellNumXInput.y + cellNumXInput.height + 1.5 * app.spacing
                width: configRectangle.width / 2.5 - app.spacing
                height: configRectangle.height / 5 - app.spacing
                validator: DoubleValidator {
                    bottom: 0.0
                    top: 10000.0
                }
                text: qsTr("0.8")
                font.family: cellNumXInput.font.family
                font.pointSize: shocks4ButtonText.font.pointSize - 5
                placeholderTextColor: "#d3d3d3"
                placeholderText: qsTr("End Time")
                color: "black"
                background: Rectangle {
                    color: "#f5f5f5"
                    radius: app.spacing / 3
                }
            }
            TextField {
                id: gammaInput
                x: cellNumYInput.x
                y: endTInput.y
                width: configRectangle.width / 2.5 - app.spacing
                height: configRectangle.height / 5 - app.spacing
                validator: DoubleValidator {
                    bottom: 1.0
                    top: 1.666665
                }
                text: qsTr("1.4")
                font.family: cellNumXInput.font.family
                font.pointSize: shocks4ButtonText.font.pointSize - 6
                placeholderTextColor: "#d3d3d3"
                placeholderText: qsTr("ɣ")
                color: "black"
                background: Rectangle {
                    color: "#f5f5f5"
                    radius: app.spacing / 3
                }
            }
            Slider {
                id: cflSlider
                from: 0
                to: 1
                value: 0.5
                stepSize: 0.01
                x: cellNumXInput.x - app.spacing / 2
                y: endTInput.y + cellNumXInput.height + app.spacing
                width: configRectangle.width - 10 * app.spacing
                Material.accent: "deepskyblue"
            }
            Rectangle {
                id: cflRectangle
                x: cflSlider.x + cflSlider.width
                y: endTInput.y + cellNumXInput.height + 1.6 * app.spacing
                width: configRectangle.width - x
                height: configRectangle.height / 20
                color: "white"
                radius: shocks4Button.radius
                Text {
                    id: cfl
                    text: qsTr("CFL: " + cflSlider.value.toFixed(2))
                    font.pointSize: shocks4ButtonText.font.pointSize - 5
                    font.family: cellNumXInput.font.family
                }
            }
            Slider {
                id: alphaSlider
                x: cflSlider.x
                y: cflSlider.y + 2.4 * app.spacing
                from: 1.0
                to: 2.0
                value: 2.0
                stepSize: 0.01
                width: configRectangle.width - 10 * app.spacing
                Material.accent: "deepskyblue"
            }
            Rectangle {
                id: alphaRectangle
                x: alphaSlider.x + alphaSlider.width + 1.1 * app.spacing
                y: endTInput.y + cellNumXInput.height + 3.9 * app.spacing
                width: configRectangle.width - x
                height: configRectangle.height / 20
                color: "white"
                radius: shocks4Button.radius
                Text {
                    id: alpha
                    text: qsTr("ɑ: " + alphaSlider.value.toFixed(2))
                    font.pointSize: shocks4ButtonText.font.pointSize - 5
                    font.family: cellNumXInput.font.family
                }
            }
            Button {
                id: buttonGO
                x: cellNumYInput.x + 0.8 * app.spacing
                y: configRectangle.height - 3.5 * app.spacing
                width: 5 * app.spacing
                contentItem: Text {
                    text: qsTr("Go")
                    font.pixelSize: shocks4ButtonText.font.pixelSize - 5
                    font.family: cellNumXInput.font.family
                    color: "black"
                    horizontalAlignment: Text.AlignHCenter
                    verticalAlignment: Text.AlignVCenter
                    anchors.fill: parent
                }
                background: Rectangle {
                    id: buttonGoRect
                    color: "#fc9b76"
                    radius: app.spacing
                }
                MouseArea {
                    anchors.fill: parent
                    hoverEnabled: true
                    onEntered: {
                        if (buttonGO.enabled)
                            buttonGoRect.color = "#ffc4a7";
                    }
                    onExited: buttonGoRect.color = "#fc9b76"
                    onClicked: KernelRelated.onGoClicked(app.testID)
                }
            }
            Button {
                id: buttonResult
                x: cellNumYInput.x - 7.5 * app.spacing
                y: configRectangle.height - 3.5 * app.spacing
                width: 5 * app.spacing
                enabled: false
                contentItem: Text {
                    text: qsTr("Result")
                    font.pixelSize: shocks4ButtonText.font.pixelSize - 5
                    font.family: cellNumXInput.font.family
                    color: "black"
                    horizontalAlignment: Text.AlignHCenter
                    verticalAlignment: Text.AlignVCenter
                    anchors.fill: parent
                }
                background: Rectangle {
                    id: buttonShowResultRect
                    color: "#1385ff"
                    radius: app.spacing
                }
                MouseArea {
                    anchors.fill: parent
                    hoverEnabled: true
                    onEntered: {
                        if (buttonResult.enabled) {
                            buttonShowResultRect.color = "#8ac1ff";
                        }
                    }
                    onExited: buttonShowResultRect.color = "#1385ff"
                    onClicked: {
                        if (buttonResult.enabled) {
                            plot.source = "";
                            plot.source = "file:./data/data.png";
                            plotWind.visible = true;
                        }
                    }
                }
            }
        }//config rectangle

        Window {
            id: plotWind
            width: mainWind.width
            height: mainWind.height
            visible: false
            Rectangle {
                anchors.fill: parent
                Image {
                    id: plot
                    anchors.fill: parent
                    anchors.margins: 5
                    source: "file:./data/data.png"
                    fillMode: Image.PreserveAspectFit
                    horizontalAlignment: Image.AlignHCenter
                    verticalAlignment: Image.AlignVCenter
                    smooth: true
                    cache: false
                }
            }
        }
    }

    Connections {
        target: launcher
        function onProgressChanged() {
            KernelRelated.onProgressChanged(app.testID, launcher.progress);
        }
        function onSimulationFinished() {
            KernelRelated.onSimulationCompleted();
        }
        function onPlotFinished(filePath) {
            KernelRelated.onPlotCompleted(filePath);
        }
    }
}//application window

