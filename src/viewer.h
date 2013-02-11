#ifndef _VIEWER_H_
#define _VIEWER_H_

// GL
#include <QGLViewer/qglviewer.h>

// local
#include "scene.h"

class Viewer : public QGLViewer{
    Q_OBJECT
    
private:
    Scene* m_scene;
    bool m_view_dual;
    bool m_view_primal;

public:
    Viewer(QWidget* parent) : QGLViewer(parent)
    {
        setBackgroundColor(::Qt::white);
        m_scene = NULL;
        m_view_dual = false;
        m_view_primal = true;
    }
    
    ~Viewer()
	{
	}
    
    Scene* get_scene() const
    {
        return m_scene;
    }
    
    void set_scene(Scene* scene)
    {
        m_scene = scene;
    }
    
    void toggle_view_primal() { m_view_primal = !m_view_primal; }
    
    void toggle_view_dual() { m_view_dual = !m_view_dual; }
    
    void initializeGL()
    {
        QGLViewer::initializeGL();
        glDisable(GL_DITHER);
        glDisable(GL_LIGHT0);
        glDisable(GL_LIGHTING);
        glDisable(GL_CULL_FACE);
        glEnable(GL_DEPTH_TEST);
        glEnable(GL_LINE_SMOOTH);
        glHint(GL_LINE_SMOOTH_HINT,GL_NICEST);
    }

    void draw()
    {
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glClearColor(1.0f,1.0f,1.0f,0.0f);
        if (!m_scene) return;
        draw_scene();
    }
    
    void draw_scene()
    {
        glPushMatrix();
        glPopMatrix();
    }
};

#endif
