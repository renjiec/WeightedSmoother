#ifndef _WINDOW_
#define _WINDOW_

// Qt
#include <QWidget>
#include <QString>

// local
#include "ui_window.h"
#include "scene.h"

class MainWindow : public QMainWindow, public Ui_MainWindow
{
	Q_OBJECT
    
private:
    Scene m_scene;
    
public:
	MainWindow();
    ~MainWindow();
    
public slots:
    
    // Drag & drop
    void dropEvent(QDropEvent *event);
    void dragEnterEvent(QDragEnterEvent *event);
    void open(QString filename);

    // View
    void on_actionViewPrimal_toggled();
    void on_actionViewDual_toggled();
    
    // File
    void on_actionOpenMesh_triggered();
    void on_actionSaveMesh_triggered();
};

#endif
