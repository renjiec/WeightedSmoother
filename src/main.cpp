#include "window.h"

#include <iostream>
#include <fstream>


int main(int argv, char **args)
{	
	QApplication app(argv, args);
	app.setApplicationName("Smoother");

	MainWindow window;
	window.show();
	return app.exec();
}
