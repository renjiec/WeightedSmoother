// Qt
#include <QtGui>
#include <QFileDialog>
#include <QInputDialog>

// local
#include "window.h"

MainWindow::MainWindow() : QMainWindow(), Ui_MainWindow()
{
    // init viewer
	setupUi(this);
	setAcceptDrops(true);
    viewer->set_scene(&m_scene);
}

MainWindow::~MainWindow()
{
}

// drwag & drop

void MainWindow::dragEnterEvent(QDragEnterEvent *event)
{
	if (event->mimeData()->hasFormat("text/uri-list"))
		event->acceptProposedAction();
}

void MainWindow::dropEvent(QDropEvent *event)
{
	Q_FOREACH(QUrl url, event->mimeData()->urls())
    {
		QString filename = url.toLocalFile();
		if(!filename.isEmpty())
		{
			QTextStream(stderr) << QString("drop event(\"%1\")\n").arg(filename);
			open(filename);
		}
	}
	event->acceptProposedAction();
}

void MainWindow::open(QString filename)
{
	QFileInfo fileinfo(filename);
	if (fileinfo.isFile() && fileinfo.isReadable())
	{
        m_scene.open_mesh(filename.toStdString());
        viewer->repaint();
	}
}

// View //

void MainWindow::on_actionViewPrimal_toggled()
{
    viewer->toggle_view_primal();
    viewer->repaint();
}

void MainWindow::on_actionViewDual_toggled()
{
    viewer->toggle_view_dual();
    viewer->repaint();
}

// File //

void MainWindow::on_actionOpenMesh_triggered()
{
    QString file = QFileDialog::getOpenFileName(this, "OFF", "../../data", "*.*");
    open(file);
}

void MainWindow::on_actionSaveMesh_triggered()
{
    QString file = QFileDialog::getSaveFileName(this, "Save", "../../data");
    if (file.isEmpty()) return;
    m_scene.save_mesh(file.toStdString());
}
